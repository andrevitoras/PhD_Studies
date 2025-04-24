import time
from datetime import datetime
from pathlib import Path

from numpy import array, zeros, arange, sqrt
from pandas import read_csv, DataFrame
from tqdm import tqdm

from solaren.niopy.geometric_transforms import dst, nrm
from solaren.pysoltrace import PySolTrace
from solaren.scopy import OpticalProperty
from solaren.scopy.linear_fresnel import (Absorber, uniform_centers, PrimaryMirror, PrimaryField,
                                          primaries_curvature_radius)

from solaren.scopy.nio_concentrators import cpc_type, symmetric_cec2tube, symmetric_cec2evacuated_tube
from solaren.scopy.sunlight import RadialSource, sun_direction, SiteData
from solaren.soltracepy import (Trace, Optics, Stage, Geometry, ElementStats, soltrace_script, run_soltrace,
                                ElementFlux, read_element_stats)
from utils import dic2json


########################################################################################################################
# CLASSES ##############################################################################################################


class OpticalSettings:

    def __init__(self,
                 rho=1.0, alpha=1.0,
                 slope_error: float = 0.,
                 specular_error: float = 0.):

        self.rho = abs(rho)
        self.alpha = abs(alpha)

        self.slope_error = 0.0 if slope_error is None else abs(slope_error)
        self.specular_error = 0.0 if specular_error is None else abs(specular_error)
        self.overall_error = sqrt(4 * self.slope_error ** 2 + self.specular_error ** 2)
        self.is_errors = False if round(self.overall_error, 10) == 0. else True

        self.primaries_property = OpticalProperty.reflector(name='primary_property',
                                                            rho=self.rho,
                                                            slope_error=self.slope_error,
                                                            spec_error=self.specular_error)

        self.secondary_property = OpticalProperty.secondary(name='secondary_property',
                                                            rho=self.rho,
                                                            slope_error=self.slope_error,
                                                            spec_error=self.specular_error)

        self.absorber_property = OpticalProperty.flat_absorber(name='absorber_property', alpha=self.alpha)
        self.tube_property = OpticalProperty.absorber_tube(name='absorber_property', alpha=self.alpha)


class SunSettings:

    def __init__(self,
                 sunshape: str = None,
                 size: float = None,
                 rays: float = 1e6,
                 cpus=10,
                 seed=123):

        # Sunshape attributes
        self.profile = sunshape
        self.sunshape_size = abs(size)

        self.source = RadialSource(profile=sunshape, size=size)
        self.is_sunshape = False if sunshape is None else True

        # Sun and Trace attributes
        self.rays = int(abs(rays))
        self.max_rays_traced = 100 * self.rays
        self.cpus = int(abs(cpus))
        self.seed = seed


class EffectiveSource:

    def __init__(self,
                 name: str,
                 sun: SunSettings, optical: OpticalSettings):
        self.name = name

        self.sun = sun
        self.optical = optical

        self.is_sunshape = self.sun.is_sunshape
        self.is_surface_errors = self.optical.is_errors

        self.full_prop = self.optical.rho**2 * self.optical.alpha
        self.prop = self.optical.rho * self.optical.alpha

        self.errors_rms_width = self.optical.overall_error * sqrt(2)
        if self.sun.is_sunshape:
            self.rms_width = sqrt(self.sun.source.rms_width**2 + self.errors_rms_width**2)
        else:
            self.rms_width = self.errors_rms_width

    def cum_eff(self):

        cum_eff = self.sun.source.linear_convolution(slope_error=self.optical.slope_error,
                                                     specular_error=self.optical.specular_error)
        return cum_eff


class uniform_geometry:

    def __init__(self,
                 name: str,
                 number_of_mirrors: int,
                 mirror_width: float,
                 mirror_shift: float,
                 receiver_height: float,
                 length: float,
                 radius_design='zenithal',
                 nbr_points=None):

        self.name = name
        self.number_primaries = abs(number_of_mirrors)
        self.mirror_width = abs(mirror_width)
        self.mirror_shift = abs(mirror_shift)
        self.receiver_height = abs(receiver_height)
        self.radius_design = radius_design

        if nbr_points is None:
            self.nbr_points = int(self.mirror_width * 300/750)
        else:
            self.nbr_points = int(nbr_points)

        self.length = abs(length)

        self.absorber = Absorber.evacuated_tube(center=array([0, receiver_height]),
                                                absorber_radius=35.0,
                                                outer_cover_radius=62.5, inner_cover_radius=62.5 - 3)

        self.centers = uniform_centers(mirror_width=self.mirror_width,
                                       center_distance=self.mirror_shift,
                                       nbr_mirrors=self.number_primaries)

        self.primary_field_width = dst(self.centers[0], self.centers[-1]) + self.mirror_width

        self.secondary = self.add_cec_secondary()

        self.rec_aim = self.secondary.aperture_center
        self.s1 = self.secondary.contour[0]
        self.s2 = self.secondary.contour[-1]

        self.flat_absorber = Absorber.flat(width=dst(self.s1, self.s2),
                                           center=self.rec_aim,
                                           axis=nrm(self.s2 - self.s1))

        self.primary_field = self.design_primary_field(radius_design=radius_design)

        self.net_area = self.number_primaries * self.mirror_width * self.length
        self.gross_area = self.primary_field_width * self.length

        self.filling_factor = self.net_area / self.gross_area

    def add_cec_secondary(self,
                          gap_radius: float = None,
                          nbr_contour_points: int = 200) -> cpc_type:
        points_per_section = int(round(nbr_contour_points / 4))
        secondary_curve = cec4tube(primary_field_width=self.primary_field_width,
                                   tube_center=self.absorber.center,
                                   tube_radius=self.absorber.absorber_tube.radius,
                                   gap_radius=self.absorber.outer_radius + 1 if gap_radius is None else gap_radius,
                                   points_per_section=points_per_section)

        return secondary_curve

    def design_primary_field(self, radius_design='zenithal'):

        radii = primaries_curvature_radius(centers=self.centers, rec_aim=self.rec_aim, radius_design=radius_design)

        heliostats = [PrimaryMirror(width=self.mirror_width, radius=r, center=hc, nbr_pts=self.nbr_points)
                      for hc, r in zip(self.centers, radii)]

        self.primary_field = PrimaryField(heliostats=heliostats)

        return self.primary_field

    def analytical_evaluation(self,
                              ES: EffectiveSource,
                              symmetric: bool,
                              path: Path,
                              force_sim: bool,
                              full_optic: bool = False,
                              file_name: str = None) -> DataFrame:

        df = analytical_analysis(lfc=self, ES=ES, symmetric=symmetric, full_optic=full_optic,
                                 path=path, force_sim=force_sim, file_name=file_name)

        return df

    def raytracing_evaluation(self,
                              ES: EffectiveSource,
                              symmetric: bool,
                              path: Path,
                              engine: str,
                              force_sim: bool,
                              full_optic: bool) -> DataFrame:

        df = raytracing_analysis(lfc=self, ES=ES, symmetric=symmetric,
                                 path=path, force_sim=force_sim,
                                 full_optic=full_optic, engine=engine)

        return df

    def validation_analysis(self,
                            ES: EffectiveSource,
                            path: Path,
                            engine: str,
                            force_sim: bool):

        # Updating the number of desired ray intersections to be traced in SolTrace.
        # The density of 2517 ray intersections per square meter of gross aperture was determined in the
        # convergence analysis
        nray_hits = int(2517 * self.gross_area / 1e6)

        # to ensure that nray_hits is never lower than 2e5 rays intersections
        nray_hits = nray_hits if nray_hits > 2e5 else 200000

        # updating the number of max traced rays
        ES.sun.max_rays_traced = 100 * ES.sun.rays

        # fixing a seed number to ensure replication of results
        ES.sun.seed = 123

        if engine == 'GUI':
            n_cpus = 10
        else:
            if nray_hits <= 1e6:
                n_cpus = 1
            elif 1e6 < nray_hits < 2.5e6:
                n_cpus = 5
            else:
                n_cpus = 10

        ES.sun.cpus = n_cpus
        ################################################################################################################

        # Analytical and ray-tracing results
        an_df = self.analytical_evaluation(ES=ES, symmetric=False, path=path,
                                           force_sim=force_sim,
                                           full_optic=False)

        rt_df = self.raytracing_evaluation(engine=engine,
                                           ES=ES, symmetric=False, path=path,
                                           force_sim=force_sim,
                                           full_optic=False)
        ################################################################################################################

        return an_df, rt_df


########################################################################################################################
# FUNCTIONS ############################################################################################################

def analytical_analysis(lfc: uniform_geometry,
                        ES: EffectiveSource,
                        symmetric: bool,
                        path: Path,
                        force_sim: bool,
                        full_optic: bool,
                        file_name=None):

    files_path = Path(path, lfc.name, ES.name)
    files_path.mkdir(parents=True, exist_ok=True)
    results_file = Path(files_path,
                        f'analytical_data.csv' if file_name is None else f'{file_name}_analytical_data.csv')

    if results_file.is_file() and not force_sim:
        df = read_csv(results_file)
    else:
        cum_eff = ES.cum_eff()
        print(f'Starting bi-axial analytical optical analysis for {lfc.name} and {ES.name} at {datetime.now()}')
        efficiency_data = lfc.primary_field.optical_analysis(flat_absorber=lfc.flat_absorber,
                                                             length=lfc.length,
                                                             rec_aim=lfc.rec_aim,
                                                             cum_eff=cum_eff,
                                                             symmetric=symmetric,
                                                             end_losses=True,
                                                             factorized=False)

        df = DataFrame(efficiency_data, columns=['theta_t', 'theta_l', 'eta'])

        df['eta'] = df['eta'] * ES.full_prop if full_optic else df['eta'] * ES.prop

        df.to_csv(results_file, index=False)

    return df


def raytracing_analysis(lfc: uniform_geometry,
                        ES: EffectiveSource,
                        path: Path,
                        full_optic: bool,
                        force_sim: bool,
                        engine: str,
                        symmetric=False):

    lv = 0. if symmetric else -90
    angles_list = array([[x, y] for x in arange(lv, 95., 5.) for y in arange(0., 95., 5.)])

    if engine in ('GUI', 'API'):
        files_path = Path(path, f'{lfc.name}_full' if full_optic else lfc.name, ES.name, engine)
        files_path.mkdir(parents=True, exist_ok=True)
    else:
        raise ValueError('Please, select and available engine to run the ray-tracing simulation')

    results_file = Path(files_path.parent, f'{engine}_raytracing_data.csv')

    if results_file.is_file() and not force_sim:
        df = read_csv(results_file)
    else:
        print(f'Starting bi-axial raytracing optical analysis for {lfc.name} and {ES.name} at {datetime.now()}')
        efficiencies = [raytracing_efficiency(lfc=lfc,
                                              theta_t=pair[0], theta_l=pair[1],
                                              ES=ES,
                                              force_sim=force_sim, full_optic=full_optic, engine=engine,
                                              files_path=files_path, file_name='optic')
                        for pair in tqdm(angles_list)]

        biaxial_data = zeros(shape=(angles_list.shape[0], 3))
        biaxial_data.T[:] = angles_list.T[0], angles_list.T[1], efficiencies

        df = DataFrame(biaxial_data, columns=['theta_t', 'theta_l', 'eta'])
        df.to_csv(results_file, index=False)

    return df


def raytracing_factorized_analysis(lfc: uniform_geometry,
                                   ES: EffectiveSource,
                                   path: Path,
                                   full_optic: bool,
                                   force_sim: bool,
                                   engine: str,
                                   symmetric=False):

    lv = 0. if symmetric else -90
    transversal_angles = array([x for x in arange(lv, 95., 5.)])
    longitudinal_angles = array([x for x in arange(0., 95., 5.)])

    if engine in ('GUI', 'API'):
        files_path = Path(path,
                          f'{lfc.name}_full' if full_optic else lfc.name,
                          ES.name,
                          engine)
        files_path.mkdir(parents=True, exist_ok=True)
    else:
        raise ValueError('Please, select and available engine to run the ray-tracing simulation')

    print(f'Starting transversal raytracing optical analysis at {datetime.now()}')
    transversal_eta = [raytracing_efficiency(lfc=lfc,
                                             theta_t=theta, theta_l=0.,
                                             ES=ES,
                                             force_sim=force_sim, full_optic=full_optic, engine=engine,
                                             files_path=files_path, file_name='optic')
                       for theta in tqdm(transversal_angles)]

    print(f'Starting transversal raytracing optical analysis at {datetime.now()}')
    longitudinal_eta = [raytracing_efficiency(lfc=lfc,
                                              theta_t=0., theta_l=theta,
                                              ES=ES,
                                              force_sim=force_sim, full_optic=full_optic, engine=engine,
                                              files_path=files_path, file_name='optic')
                        for theta in tqdm(longitudinal_angles)]

    transversal_data = zeros(shape=(transversal_angles.shape[0], 2))
    longitudinal_data = zeros(shape=(longitudinal_angles.shape[0], 2))

    transversal_data.T[:] = transversal_angles, transversal_eta
    longitudinal_data.T[:] = longitudinal_angles, longitudinal_eta

    return transversal_data, longitudinal_data


def raytracing_efficiency(lfc: uniform_geometry,
                          theta_t: float, theta_l: float,
                          ES: EffectiveSource,
                          full_optic: bool,
                          files_path: Path, file_name: str,
                          force_sim: bool, engine: str):

    if engine == 'GUI':
        optical_efficiency = lfc_gui_efficiency(lfc=lfc,
                                                theta_t=theta_t, theta_l=theta_l,
                                                ES=ES,
                                                full_optic=full_optic,
                                                file_name=file_name, files_path=files_path,
                                                force_sim=force_sim)
    elif engine == 'API':
        optical_efficiency = lfc_api_efficiency(lfc=lfc,
                                                theta_t=theta_t, theta_l=theta_l,
                                                ES=ES,
                                                full_optic=full_optic,
                                                file_name=file_name, files_path=files_path,
                                                force_sim=force_sim)
    else:
        raise ValueError('Please, select and available engine to run the ray-tracing simulation')

    return optical_efficiency


def lfc_gui_efficiency(lfc: uniform_geometry,
                       theta_t: float, theta_l: float,
                       ES: EffectiveSource,
                       full_optic: bool,
                       files_path: Path, file_name: str, force_sim: bool):

    # By definition, if any of the incidence angle is 90º, efficiency is set to zero.
    if abs(theta_t) == 90. or abs(theta_l) == 90.:
        optical_efficiency = 0.
    # If angles are not 90º, run calculations
    else:

        # Incidence sunlight direction
        sun_dir = sun_direction(theta_t=theta_t, theta_l=theta_l)

        # Optics
        if full_optic:
            optics = [ES.optical.primaries_property.to_soltrace(),
                      ES.optical.secondary_property.to_soltrace(),
                      ES.optical.tube_property.to_soltrace()]
        else:
            optics = [ES.optical.primaries_property.to_soltrace(), ES.optical.absorber_property.to_soltrace()]

        # Elements
        # primary field
        elements = lfc.primary_field.to_soltrace(rec_aim=lfc.rec_aim, sun_dir=sun_dir,
                                                 optic=optics[0], length=lfc.length)

        # Secondary and absorber tube
        if full_optic:
            elements += lfc.secondary.as_plane_curve().to_soltrace(name='secondary_optic', length=lfc.length,
                                                                   optic=optics[1])
            elements += lfc.absorber.absorber_tube.as_soltrace_element(length=lfc.length,
                                                                       optic=optics[2], name='absorber_tube')
        # Secondary aperture as a flat absorber
        else:
            elements += lfc.flat_absorber.as_soltrace_element(length=lfc.length, optic=optics[-1],
                                                              name='flat_absorber')

        # Creating SolTrace objects for main input Boxes
        sun = ES.sun.source.to_soltrace(sun_dir=sun_dir)
        optics = Optics(properties=optics)
        stage = Stage(name='linear_fresnel', elements=elements)
        geometry = Geometry(stages=[stage])

        stats = ElementStats(stats_name=f'{file_name}_absorber_flux_{theta_t}_{theta_l}',
                             stage_index=0, element_index=len(elements) - 1,
                             dni=1000., x_bins=5, y_bins=5, final_rays=True)

        trace_options = Trace(rays=ES.sun.rays, cpus=ES.sun.cpus, seed=ES.sun.seed,
                              sunshape=ES.is_sunshape, optical_errors=ES.is_surface_errors,
                              point_focus=False, simulate=True)

        # Creating and writing the codes in a LK script file #########################################
        script_full_path = soltrace_script(file_path=files_path,
                                           file_name=f'{file_name}_{theta_t}_{theta_l}',
                                           sun=sun, optics=optics,
                                           geometry=geometry,
                                           trace=trace_options, stats=stats)

        # Running the script file with the SolTrace 2012.7.9 version from the prompt #################
        if not stats.file_full_path.is_file() or force_sim:
            t0 = time.time()
            run_soltrace(script_full_path)
            sim_time = time.time() - t0

            # Reading exported file by SolTrace and add the simulation run time to the dict, and then exporting it
            absorber_stats = read_element_stats(stats.file_full_path)
            absorber_stats['sim_time'] = sim_time
            dic2json(d=absorber_stats,
                     file_name=stats.file_full_path.name[:-5],
                     file_path=Path(stats.file_full_path.parent), indent=4)

            absorber_stats = ElementFlux(stats.file_full_path)
        else:
            absorber_stats = ElementFlux(stats.file_full_path)

        # Reading the stats json file with the main absorber flux data #################################################
        # It then calculates and returns the optical efficiency
        # Here, optical efficiency is defined as the ratio between the absorber flux at the receiver and the
        # incident flux at the mirrors aperture as if the sun was perpendicular to each mirror.

        # OBS: 'linear_fresnel' module measures are in millimeters, and SolTrace works with meters.
        # absorber_stats = ElementFlux(stats.file_full_path)
        # absorber_flux = absorber_stats.total_flux * 1000.

        absorber_flux = absorber_stats.total_flux * 1000.
        optical_efficiency = absorber_flux / (stats.dni * lfc.primary_field.widths.sum() * lfc.length * 1e-6)

        ################################################################################################################

    return optical_efficiency


def lfc_api_efficiency(lfc: uniform_geometry,
                       theta_t: float, theta_l: float,
                       ES: EffectiveSource,
                       full_optic: bool,
                       files_path: Path, file_name: str, force_sim: bool):

    # By definition, if any of the incidence angle is 90º, efficiency is set to zero.
    if abs(theta_t) == 90. or abs(theta_l) == 90.:
        optical_efficiency = 0.
    # If angles are not 90º, run calculations
    else:

        # The PySolTrace object ############
        optical_system = PySolTrace()

        ####################################

        # Sun ###############################################################
        sun_dir = sun_direction(theta_t=theta_t,
                                theta_l=theta_l)
        optical_system.sun = ES.sun.source.to_pysoltrace(sun_dir=sun_dir)
        #####################################################################

        # Optics ########################################################################
        if full_optic:
            optics = [ES.optical.primaries_property.to_pysoltrace(id_number=0),
                      ES.optical.secondary_property.to_pysoltrace(id_number=1),
                      ES.optical.tube_property.to_pysoltrace(id_number=2)]
        else:
            optics = [ES.optical.primaries_property.to_pysoltrace(id_number=0),
                      ES.optical.absorber_property.to_pysoltrace(id_number=1)]

        # updating the optics of the PySolTrace object -- Optics Box
        optical_system.optics = optics
        ################################################################################

        # Stage #########################################################################################
        lfc_stage = PySolTrace.Stage(id=0)
        lfc_stage.name = 'linear_fresnel'

        # Elements
        # primary field
        elements = [hel.to_pysoltrace(sun_dir=sun_dir,
                                      rec_aim=lfc.rec_aim,
                                      length=lfc.length,
                                      parent_stage=lfc_stage,
                                      id_number=i,
                                      optic=optics[0])
                    for i, hel in enumerate(lfc.primary_field.heliostats)]

        # Secondary optic and absorber tube
        if full_optic:
            elements += lfc.secondary.as_plane_curve().to_pysoltrace(length=lfc.length,
                                                                     parent_stage=lfc_stage,
                                                                     id_number=len(elements),
                                                                     optic=optics[1])

            elements += [lfc.absorber.absorber_tube.to_pysoltrace(length=lfc.length,
                                                                  parent_stage=lfc_stage,
                                                                  id_number=len(elements),
                                                                  optic=optics[2])]
        # Secondary aperture as a flat absorber
        else:
            elements += [lfc.flat_absorber.to_pysoltrace(length=lfc.length,
                                                         parent_stage=lfc_stage,
                                                         id_number=len(elements),
                                                         optic=optics[1])]
        # updating the elements in the stage
        lfc_stage.elements = elements

        # updating the stages of the PySolTrace object -- Geometry Box
        optical_system.stages = [lfc_stage]

        # Ray-tracing settings
        optical_system.num_ray_hits = ES.sun.rays
        optical_system.max_rays_traced = ES.sun.max_rays_traced

        optical_system.is_sunshape = ES.is_sunshape
        optical_system.is_surface_errors = ES.is_surface_errors

        #################################################################################################

        input_file = Path(files_path, f'{file_name}_{theta_t}_{theta_l}.stinput')
        stats_file = Path(files_path, f'{file_name}_{theta_t}_{theta_l}_absorber_flux_stats.json')
        if stats_file.is_file() and not force_sim:
            absorber_stats = ElementFlux(stats_file)

        else:
            # absorber stats file does not exist. Run simulation and export the dict as json file
            t0 = time.time()
            optical_system.run(seed=ES.sun.seed, nthread=ES.sun.cpus, as_power_tower=False)
            optical_system.write_soltrace_input_file(path=str(input_file))
            absorber_stats = optical_system.element_stats(element=optical_system.stages[0].elements[-1],
                                                          xbins=5, ybins=5,
                                                          absorbed_only=True)
            absorber_stats['sim_time'] = time.time() - t0
            dic2json(d=absorber_stats,
                     file_path=files_path,
                     file_name=f'{file_name}_{theta_t}_{theta_l}_absorber_flux_stats', indent=4)
            absorber_stats = ElementFlux(stats_file)

        absorber_flux = absorber_stats.total_flux * 1000.
        optical_efficiency = absorber_flux / (optical_system.dni * lfc.primary_field.widths.sum() * lfc.length * 1e-6)

    return optical_efficiency


def cec4tube(primary_field_width: float,
             tube_center: array,
             tube_radius: float,
             gap_radius: float = None,
             points_per_section: int = 50) -> cpc_type:
    p_left_edge, p_right_edge = array([-primary_field_width / 2, 0]), array([primary_field_width / 2, 0])
    source_width = dst(p_left_edge, p_right_edge)
    source_distance = (tube_center - p_left_edge)[-1]

    # The symmetric CEC secondary with no gap between the optic and absorber tube.
    if abs(gap_radius) == 0 or gap_radius is None:
        l_con, l_inv, r_inv, r_con = symmetric_cec2tube(tube_radius=abs(tube_radius),
                                                        tube_center=tube_center,
                                                        source_width=source_width,
                                                        source_distance=source_distance,
                                                        nbr_pts=points_per_section, upwards=False)
    # The symmetric CEC secondary with a gap between the optic and absorber tube.
    elif abs(gap_radius) > abs(tube_radius):
        l_con, l_inv, r_inv, r_con = symmetric_cec2evacuated_tube(tube_radius=abs(tube_radius),
                                                                  tube_center=tube_center,
                                                                  cover_radius=gap_radius, dy=0,
                                                                  source_width=source_width,
                                                                  source_distance=source_distance,
                                                                  nbr_pts=points_per_section, upwards=False)
    else:
        raise ValueError('The gap and tube radii do not match: tube radius should be lower than gap radius!')

    optic = cpc_type(left_conic=l_con, left_involute=l_inv, right_involute=r_inv, right_conic=r_con)

    return optic
