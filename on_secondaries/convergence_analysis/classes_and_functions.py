import datetime
from numpy import rad2deg, zeros

from solaren.niopy.geometric_transforms import dst, tgs2tube, ang
from solaren.scopy import OpticalProperty
from solaren.scopy.linear_fresnel import (PrimaryField, Absorber, uniform_centers,
                                          primaries_curvature_radius, PrimaryMirror)
from solaren.scopy.nio_concentrators import symmetric_cec2evacuated_tube, cpc_type, symmetric_cpc2evacuated_tube
from solaren.scopy.sunlight import RadialSource, sun_direction
from solaren.soltracepy import *


class OpticalSettings:

    def __init__(self,
                 primaries_property: OpticalProperty.reflector,
                 secondary_property: OpticalProperty.secondary,
                 absorber_property: OpticalProperty.evacuated_tube):

        self.primaries_property = primaries_property
        self.secondary_property = secondary_property
        self.absorber_property = absorber_property

        self.primaries_prop = primaries_property.to_soltrace()
        self.secondaries_prop = secondary_property.to_soltrace()
        self.o_cover_prop, self.i_cover_prop, self.abs_prop = self.absorber_property.to_soltrace()

        self.properties = [self.primaries_prop,
                           self.secondaries_prop,
                           self.o_cover_prop, self.i_cover_prop, self.abs_prop]


class lfc_optic:

    def __init__(self,
                 name: str,
                 mirror_width: float, center_distance: float,
                 nbr_mirrors: int,
                 receiver_height: float,
                 absorber_radius: float, outer_cover_radius: float, inner_cover_radius: float,
                 concentrator_length=8000,
                 secondary_type='CEC'):

        # Main attributes ##################################################################
        self.name = name
        self.mirror_width = abs(mirror_width)
        self.secondary_type = secondary_type
        self.concentrator_length = abs(concentrator_length)

        # primary mirrors aim point is first set to 'None'.
        # it is later defined as with the secondary optic -- its aperture center point
        self.rec_aim = None
        ####################################################################################

        # The evacuated absorber tube #####################################################
        self.absorber = Absorber.evacuated_tube(center=array([0, receiver_height]),
                                                absorber_radius=absorber_radius,
                                                outer_cover_radius=outer_cover_radius,
                                                inner_cover_radius=inner_cover_radius)
        ####################################################################################

        # The center points of the primary field #################################
        self.primary_centers = uniform_centers(mirror_width=mirror_width,
                                               center_distance=center_distance,
                                               nbr_mirrors=nbr_mirrors)
        ##########################################################################

        # Primary field overall measure #############################################
        self.total_primary_width = dst(self.primary_centers[0],
                                       self.primary_centers[-1]) + mirror_width

        self.f1 = array([-0.5 * self.total_primary_width, 0])
        self.f2 = array([+0.5 * self.total_primary_width, 0])
        #############################################################################

        # The secondary optic ##############################################################
        # secondary optic design... attribute 'secondary' is a PlaneCurve object
        if secondary_type == 'CEC':
            self.secondary = self.cec_secondary()
        elif secondary_type == 'CPC':
            self.secondary = self.cpc_secondary()
        else:
            raise ValueError('Input a valid argument of secondary optic type!')
        ####################################################################################

        # Design of the primary field #######################################################

        # curvature radius of the primary mirrors
        # it considers a zenithal reference
        self.radii = primaries_curvature_radius(centers=self.primary_centers,
                                                rec_aim=self.rec_aim,
                                                radius_design='zenithal')

        # primary mirrors as 'Heliostat' objects
        heliostats = [PrimaryMirror(center=hc, width=self.mirror_width, radius=r)
                      for hc, r in zip(self.primary_centers, self.radii)
                      ]

        # primary field as 'PrimaryField' object
        self.primary_field = PrimaryField(heliostats=heliostats)
        ####################################################################################

    def tube_edges(self):

        # Right edge of the receiver  ##############################################
        # it considers the left edge of the primary field
        ta, tb = tgs2tube(point=self.f1,
                          tube_radius=self.absorber.absorber_tube.radius,
                          tube_center=self.absorber.center)

        if ta[0] < tb[0]:
            t2 = tb
        else:
            t2 = ta
        ############################################################################

        # Left edge of the receiver #################################################
        # it considers the right edge of the primary field

        ta, tb = tgs2tube(point=self.f2,
                          tube_radius=self.absorber.absorber_tube.radius,
                          tube_center=self.absorber.center)

        if ta[0] < tb[0]:
            t1 = ta
        else:
            t1 = tb
        ############################################################################

        return t1, t2

    def edge_rays_acceptance(self) -> float:

        t1, t2 = tgs2tube(point=self.f1,
                          tube_radius=self.absorber.absorber_tube.radius,
                          tube_center=self.absorber.center)

        if t1[0] < t2[0]:
            t1, t2 = t2, t1

        half_acceptance = rad2deg(ang(t1 - self.f1, array([0, 1])))

        return half_acceptance

    def cec_secondary(self, additional_gap=0., points_per_section: int = 30):

        optic_data = symmetric_cec2evacuated_tube(tube_center=self.absorber.absorber_tube.center,
                                                  tube_radius=self.absorber.absorber_tube.radius,
                                                  cover_radius=self.absorber.outer_radius,
                                                  source_distance=self.absorber.center[1],
                                                  source_width=self.total_primary_width,
                                                  nbr_pts=points_per_section, upwards=False, dy=additional_gap)

        secondary_curve = cpc_type(left_conic=optic_data[0],
                                   left_involute=optic_data[1],
                                   right_involute=optic_data[2],
                                   right_conic=optic_data[3])

        self.rec_aim = secondary_curve.aperture_center

        return secondary_curve.as_plane_curve()

    def cpc_secondary(self, acceptance_angle: float = None, additional_gap=0., points_per_section: int = 50):

        theta_a = self.edge_rays_acceptance() if acceptance_angle is None else acceptance_angle

        optic_data = symmetric_cpc2evacuated_tube(theta_a=theta_a,
                                                  tube_center=self.absorber.absorber_tube.center,
                                                  tube_radius=self.absorber.absorber_tube.radius,
                                                  cover_radius=self.absorber.outer_radius,
                                                  nbr_pts=points_per_section, upwards=False, dy=additional_gap)

        secondary_curve = cpc_type(left_conic=optic_data[0],
                                   left_involute=optic_data[1],
                                   right_involute=optic_data[2],
                                   right_conic=optic_data[3])

        self.rec_aim = secondary_curve.aperture_center

        return secondary_curve.as_plane_curve()

    def to_soltrace(self, sun_dir: array, optics: OpticalSettings) -> list:

        # The primary mirrors as SolTrace elements
        elements = self.primary_field.to_soltrace(rec_aim=self.rec_aim,
                                                  sun_dir=sun_dir,
                                                  length=self.concentrator_length,
                                                  optic=optics.primaries_prop)

        # The secondary optic as SolTrace elements
        elements += self.secondary.to_soltrace(name='secondary_optic',
                                               length=self.concentrator_length,
                                               optic=optics.secondaries_prop)

        # the evacuated absorber as SolTrace elements
        elements += self.absorber.as_soltrace_element(length=self.concentrator_length,
                                                      outer_cover_optic=optics.o_cover_prop,
                                                      inner_cover_optic=optics.i_cover_prop,
                                                      absorber_optic=optics.abs_prop)

        return elements


def simulate_optic(lfc: lfc_optic,
                   theta_t: float, theta_l: float, dni: float,
                   file_name: str, file_path: Path,
                   x_bins: int, y_bins: int,
                   source: RadialSource,
                   optics: OpticalSettings,
                   trace: Trace,
                   soltrace_version=2012,
                   force_sim=False) -> ElementFlux:

    # Ensuring that argument 'file_path' exist #######
    files_path = Path(file_path, lfc.name)
    files_path.mkdir(parents=True, exist_ok=True)
    ##################################################

    # The Sun #########################################################
    # the sun direction as a [x,y,z] vector
    sun_dir = sun_direction(theta_t=theta_t, theta_l=theta_l)
    # the object to represent Soltrace Sun box
    sun_box = source.to_soltrace(sun_dir=sun_dir)
    ###################################################################

    # The Optics
    optics_box = Optics(properties=optics.properties)

    # The Geometry
    elements = lfc.to_soltrace(sun_dir=sun_dir, optics=optics)
    stage = Stage(name='linear_fresnel', elements=elements)
    geometry_box = Geometry(stages=[stage])

    # The Element intersections stats
    absorber_stats = ElementStats(stats_name=f'absorber_flux_{file_name}_{theta_t}_{theta_l}',
                                  stage_index=0,
                                  element_index=len(elements) - 1,
                                  dni=dni,
                                  x_bins=x_bins,
                                  y_bins=y_bins,
                                  final_rays=True)

    script_path = soltrace_script(sun=sun_box,
                                  optics=optics_box,
                                  geometry=geometry_box,
                                  trace=trace,
                                  stats=absorber_stats,
                                  file_name=f'{file_name}_{theta_t}_{theta_l}',
                                  file_path=files_path,
                                  version=soltrace_version)

    # Run SolTrace from the windows prompt cmd ####################################################
    # Only the 2012.7.9 version of SolTrace can run a lk script from the prompt
    # If the absorber_stats file already exist, it should not run SolTrace
    # Finally, if forced, it will always run the simulation
    if (soltrace_version == 2012 and not absorber_stats.file_full_path.is_file()) or force_sim:
        run_soltrace(lk_file_full_path=script_path)
    ################################################################################################

    absorber_flux_data = ElementFlux(stats_file=absorber_stats.file_full_path)

    return absorber_flux_data


def get_efficiency(lfc: lfc_optic,
                   theta_t: float, theta_l: float, dni: float,
                   file_name: str, file_path: Path,
                   x_bins: int, y_bins: int,
                   source: RadialSource,
                   optics: OpticalSettings,
                   trace: Trace):

    absorber_flux_data = simulate_optic(file_name=file_name, file_path=file_path,
                                        theta_t=theta_t, theta_l=theta_l, dni=dni,
                                        lfc=lfc,
                                        x_bins=x_bins, y_bins=y_bins,
                                        source=source, optics=optics, trace=trace,
                                        soltrace_version=3)

    mirror_area = (lfc.primary_field.widths.sum() * lfc.concentrator_length) / 1e6  # im m2
    absorbed_flux = absorber_flux_data.total_flux * 1e3  # in W

    return absorbed_flux / (mirror_area * dni)


def get_uniformity(lfc: lfc_optic,
                   theta_t: float, theta_l: float, dni: float,
                   file_name: str, file_path: Path,
                   x_bins: int, y_bins: int,
                   source: RadialSource,
                   optics: OpticalSettings,
                   trace: Trace):

    absorber_flux_data = simulate_optic(file_name=file_name, file_path=file_path,
                                        theta_t=theta_t, theta_l=theta_l, dni=dni,
                                        lfc=lfc,
                                        x_bins=x_bins, y_bins=y_bins,
                                        source=source, optics=optics, trace=trace,
                                        soltrace_version=3)

    return absorber_flux_data.uniformity


def get_circumferential_uniformity(lfc: lfc_optic,
                                   theta_t: float, theta_l: float, dni: float,
                                   file_name: str, file_path: Path,
                                   x_bins: int, y_bins: int,
                                   source: RadialSource,
                                   optics: OpticalSettings,
                                   trace: Trace):

    absorber_flux_data = simulate_optic(file_name=file_name, file_path=file_path,
                                        theta_t=theta_t, theta_l=theta_l, dni=dni,
                                        lfc=lfc,
                                        x_bins=x_bins, y_bins=y_bins,
                                        source=source, optics=optics, trace=trace,
                                        soltrace_version=3)

    return absorber_flux_data.flux_x_uniformity


def get_longitudinal_uniformity(lfc: lfc_optic,
                                theta_t: float, theta_l: float, dni: float,
                                file_name: str, file_path: Path,
                                x_bins: int, y_bins: int,
                                source: RadialSource,
                                optics: OpticalSettings,
                                trace: Trace):

    absorber_flux_data = simulate_optic(file_name=file_name, file_path=file_path,
                                        theta_t=theta_t, theta_l=theta_l, dni=dni,
                                        lfc=lfc,
                                        x_bins=x_bins, y_bins=y_bins,
                                        source=source, optics=optics, trace=trace,
                                        soltrace_version=3)

    return absorber_flux_data.flux_y_uniformity


def sensitivity_analysis(rays_to_simulate: array, simulations_per_ray: int,
                         theta_t: float, dni: float,
                         lfc: lfc_optic,
                         optics: OpticalSettings,
                         source: RadialSource,
                         file_name: str,
                         file_path: Path):

    efficiency_results = zeros(shape=(rays_to_simulate.shape[0], simulations_per_ray))
    uniformity_results = zeros(shape=(rays_to_simulate.shape[0], simulations_per_ray))

    circumferential_uniformity = zeros(shape=(rays_to_simulate.shape[0], simulations_per_ray))
    longitudinal_uniformity = zeros(shape=(rays_to_simulate.shape[0], simulations_per_ray))

    files_path = Path(file_path, 'sensitivity_analysis')
    files_path.mkdir(parents=True, exist_ok=True)

    # A function to create the trace options given a particular number of rays
    def trace_options(rays: int):
        return Trace(rays=rays, seed=-1, cpus=8,
                     sunshape=True, optical_errors=True, point_focus=False, simulate=True)

    for i, nbr_rays in enumerate(rays_to_simulate):
        now = datetime.datetime.now()
        print(f'It is now {now}, and simulations for {int(nbr_rays)} rays intersections began.')

        efficiency_results[i] = [get_efficiency(theta_t=theta_t, theta_l=0., dni=dni,
                                                lfc=lfc,
                                                x_bins=150, y_bins=150,
                                                source=source, optics=optics,
                                                trace=trace_options(rays=nbr_rays),
                                                file_name=f'{file_name}_{nbr_rays}_{j + 1}',
                                                file_path=files_path) for j in range(simulations_per_ray)]

        uniformity_results[i] = [get_uniformity(theta_t=theta_t, theta_l=0., dni=dni,
                                                lfc=lfc,
                                                x_bins=150, y_bins=150,
                                                source=source, optics=optics,
                                                trace=trace_options(rays=nbr_rays),
                                                file_name=f'{file_name}_{nbr_rays}_{j + 1}',
                                                file_path=files_path) for j in range(simulations_per_ray)]

        circumferential_uniformity[i] = [get_circumferential_uniformity(theta_t=theta_t, theta_l=0., dni=dni,
                                                                        lfc=lfc,
                                                                        x_bins=150, y_bins=150,
                                                                        source=source, optics=optics,
                                                                        trace=trace_options(rays=nbr_rays),
                                                                        file_name=f'{file_name}_{nbr_rays}_{j + 1}',
                                                                        file_path=files_path)
                                         for j in range(simulations_per_ray)]

        longitudinal_uniformity[i] = [get_longitudinal_uniformity(theta_t=theta_t, theta_l=0., dni=dni,
                                                                  lfc=lfc,
                                                                  x_bins=150, y_bins=150,
                                                                  source=source, optics=optics,
                                                                  trace=trace_options(rays=nbr_rays),
                                                                  file_name=f'{file_name}_{nbr_rays}_{j + 1}',
                                                                  file_path=files_path)
                                      for j in range(simulations_per_ray)]

    return efficiency_results, uniformity_results, circumferential_uniformity, longitudinal_uniformity


def generate_script(file_name: str, file_path: Path,
                    theta_t: float, theta_l: float, dni: float,
                    lfc: lfc_optic, x_bins: int, y_bins: int,
                    source: RadialSource,
                    optics: OpticalSettings,
                    trace: Trace,
                    soltrace_version=3.0):

    # Ensuring that argument 'file_path' exist #######
    files_path = Path(file_path, lfc.name)
    files_path.mkdir(parents=True, exist_ok=True)
    ##################################################

    # The Sun #########################################################
    # the sun direction as a [x,y,z] vector
    sun_dir = sun_direction(theta_t=theta_t, theta_l=theta_l)
    # the object to represent Soltrace Sun box
    sun_box = source.to_soltrace(sun_dir=sun_dir)
    ###################################################################

    # The Optics
    optics_box = Optics(properties=optics.properties)

    # The Geometry
    elements = lfc.to_soltrace(sun_dir=sun_dir, optics=optics)
    stage = Stage(name='linear_fresnel', elements=elements)
    geometry_box = Geometry(stages=[stage])

    # The Element intersections stats
    absorber_stats = ElementStats(stats_name=f'absorber_flux_{file_name}_{theta_t}_{theta_l}',
                                  stage_index=0,
                                  element_index=len(elements) - 1,
                                  dni=dni,
                                  x_bins=x_bins,
                                  y_bins=y_bins,
                                  final_rays=True)

    script_path = soltrace_script(sun=sun_box,
                                  optics=optics_box,
                                  geometry=geometry_box,
                                  trace=trace,
                                  stats=absorber_stats,
                                  file_name=f'{file_name}_{theta_t}_{theta_l}',
                                  file_path=files_path,
                                  version=soltrace_version)

    return script_path


def get_script_code(script: Path):
    with open(script, 'r') as source_file:
        contents = source_file.readlines()

    return contents


def concatenate_scripts(scripts_list: list,
                        file_name: str, file_path: Path):

    main_file = Path(file_path, f'{file_name}.lk')

    with open(main_file, 'w') as dest_file:
        for source_file in scripts_list:

            contents = get_script_code(script=source_file)
            dest_file.writelines(contents)

    return main_file


def script4rays(nbr_rays: int, nbr_simulations: int,
                file_name: str, file_path: Path,
                theta_t: float, theta_l: float, dni: float,
                lfc: lfc_optic, x_bins: int, y_bins: int,
                source: RadialSource,
                optics: OpticalSettings,
                soltrace_version=3.0):

    files_path = Path(file_path, 'sensitivity_analysis')
    file_path.mkdir(parents=True, exist_ok=True)

    # A function to create the trace options given a particular number of rays
    def trace_options(rays: int):
        return Trace(rays=rays, seed=-1, cpus=8,
                     sunshape=True, optical_errors=True, point_focus=False, simulate=True)

    script_list = [generate_script(file_name=f'{file_name}_{nbr_rays}_{i + 1}', file_path=files_path,
                                   theta_t=theta_t, theta_l=theta_l, dni=dni,
                                   lfc=lfc, x_bins=x_bins, y_bins=y_bins,
                                   source=source, optics=optics, trace=trace_options(rays=nbr_rays),
                                   soltrace_version=soltrace_version)

                   for i in range(nbr_simulations)]

    main_script = concatenate_scripts(scripts_list=script_list,
                                      file_name=f'{file_name}_{nbr_rays}', file_path=files_path)

    return main_script


def script4sensitivity(rays_to_simulate: array, nbr_simulations: int,
                       file_name: str, file_path: Path,
                       theta_t: float, theta_l: float, dni: float,
                       lfc: lfc_optic, x_bins: int, y_bins: int,
                       source: RadialSource,
                       optics: OpticalSettings,
                       soltrace_version=3.0):

    scripts_list = [script4rays(nbr_rays=r, nbr_simulations=nbr_simulations,
                                file_name=file_name, file_path=file_path,
                                theta_t=theta_t, theta_l=theta_l, dni=dni,
                                lfc=lfc, x_bins=x_bins, y_bins=y_bins,
                                source=source, optics=optics, soltrace_version=soltrace_version)

                    for r in rays_to_simulate]

    full_script = concatenate_scripts(scripts_list=scripts_list,
                                      file_name=f'{file_name}_sensitivity_analysis', file_path=file_path)

    return full_script







