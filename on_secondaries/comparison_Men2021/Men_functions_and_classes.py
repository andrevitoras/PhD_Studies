import time
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Tuple

from numpy import array, ones, arange, zeros, absolute
from scipy.interpolate import interp1d
from pandas import read_csv
from tqdm import tqdm

from solaren.niopy.geometric_transforms import dst, tgs2tube, mid_point
from solaren.niopy.plane_curves import PlaneCurve
from solaren.scopy import OpticalProperty
from solaren.scopy.linear_fresnel import Absorber, uniform_centers, PrimaryMirror, PrimaryField, acceptance_angle
from solaren.scopy.linear_fresnel.secondaries import oommen_cpc4tube, half_acceptance2tube, cpc4tube, cec4tube
from solaren.scopy.sunlight import RadialSource, sun_direction, SiteData
from solaren.soltracepy import Trace, Optics, Stage, Geometry, ElementStats, soltrace_script, run_soltrace, ElementFlux

path_to_this_file = Path(__file__).parents
men_database_path = Path(path_to_this_file[0], f'literature_data\Men2021_raw_data.csv')


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
                 nbr_mirrors: int,
                 mirror_width: float, center_distance: float, mirror_radius: float,
                 receiver_height: float,
                 absorber_radius: float, outer_cover_radius: float, inner_cover_radius: float,
                 concentrator_length=6000):

        # Main attributes ##################################################################
        self.name = name
        self.mirror_width = abs(mirror_width)
        self.mirror_radius = abs(mirror_radius)
        self.concentrator_length = abs(concentrator_length)

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

        # Design of the primary field #######################################################
        # curvature radius of the primary mirrors
        # it considers a zenithal reference
        self.radii = ones(self.primary_centers.shape[0]) * self.mirror_radius

        # primary mirrors as 'Heliostat' objects
        heliostats = [PrimaryMirror(center=hc, width=self.mirror_width, radius=r)
                      for hc, r in zip(self.primary_centers, self.radii)]

        # primary field as 'PrimaryField' object
        self.primary_field = PrimaryField(heliostats=heliostats)
        ####################################################################################

        self.t1, self.t2 = self.tube_edges()

        self.rec_aim = self.absorber.center
        self.secondary_type = None
        self.secondary = None

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

    def add_oommen_cpc(self, theta_a: float, theta_max: float, gap_radius: float, nbr_contour_points: int = 200):

        pts_per_side = int(round(nbr_contour_points / 2))

        cpc_curve = oommen_cpc4tube(tube_center=self.absorber.center,
                                    tube_radius=self.absorber.absorber_tube.radius,
                                    gap_radius=gap_radius,
                                    theta_a=theta_a, theta_max=theta_max,
                                    points_per_side=pts_per_side)

        self.rec_aim = mid_point(p=cpc_curve.curve_pts[0], q=cpc_curve.curve_pts[-1])
        self.secondary = cpc_curve

        return cpc_curve

    def edge_rays_acceptance(self) -> float:

        half_acceptance = half_acceptance2tube(primary_field=self.primary_field.primaries,
                                               tube_center=self.absorber.center,
                                               tube_radius=self.absorber.absorber_tube.radius)

        return half_acceptance

    def add_cpc_secondary(self, gap_radius: float, nbr_contour_points: int = 200):

        points_per_section = int(round(nbr_contour_points / 4))
        secondary_curve = cpc4tube(primary_field=self.primary_field.primaries,
                                   tube_center=self.absorber.center,
                                   tube_radius=self.absorber.absorber_tube.radius,
                                   gap_radius=gap_radius,
                                   points_per_section=points_per_section)

        self.rec_aim = secondary_curve.aperture_center
        self.secondary = secondary_curve.as_plane_curve()

        return self.secondary

    def add_cec_secondary(self, gap_radius: float, nbr_contour_points: int = 200):

        points_per_section = int(round(nbr_contour_points / 4))
        secondary_curve = cec4tube(primary_field=self.primary_field.primaries,
                                   tube_center=self.absorber.center,
                                   tube_radius=self.absorber.absorber_tube.radius,
                                   gap_radius=gap_radius,
                                   points_per_section=points_per_section)

        self.rec_aim = secondary_curve.aperture_center
        self.secondary = secondary_curve.as_plane_curve()

        return self.secondary

    def to_soltrace(self, sun_dir: array, optics: OpticalSettings,
                    rec_dx=0., rec_dy=0.) -> list:

        # The primary mirrors as SolTrace elements
        elements = self.primary_field.to_soltrace(rec_aim=self.rec_aim,
                                                  sun_dir=sun_dir,
                                                  length=self.concentrator_length,
                                                  optic=optics.primaries_prop)

        # The receiver has no position error
        # Redundant codes, but it keeps the original one
        if rec_dx == 0. and rec_dy == 0.:
            # The secondary optic as SolTrace elements
            if self.secondary is not None:
                elements += self.secondary.to_soltrace(name='secondary_optic',
                                                       length=self.concentrator_length,
                                                       optic=optics.secondaries_prop)
            else:
                raise ValueError('No secondary optic was identified. Please add one!')

            # the evacuated absorber as SolTrace elements
            elements += self.absorber.as_soltrace_element(length=self.concentrator_length,
                                                          outer_cover_optic=optics.o_cover_prop,
                                                          inner_cover_optic=optics.i_cover_prop,
                                                          absorber_optic=optics.abs_prop)
        # The receiver has a position error
        else:
            if self.secondary is not None:
                displaced_secondary = PlaneCurve(curve_pts=self.secondary.curve_pts + array([rec_dx, rec_dy]),
                                                 curve_center=self.secondary.center + array([rec_dx, rec_dy]))

                elements += displaced_secondary.to_soltrace(name='secondary_optic',
                                                            length=self.concentrator_length,
                                                            optic=optics.secondaries_prop)
            else:
                raise ValueError('No secondary optic was identified. Please add one!')

            displaced_absorber = Absorber.evacuated_tube(center=self.absorber.center + array([rec_dx, rec_dy]),
                                                         absorber_radius=self.absorber.radius,
                                                         outer_cover_radius=self.absorber.outer_radius,
                                                         inner_cover_radius=self.absorber.inner_radius)

            elements += displaced_absorber.as_soltrace_element(length=self.concentrator_length,
                                                               outer_cover_optic=optics.o_cover_prop,
                                                               inner_cover_optic=optics.i_cover_prop,
                                                               absorber_optic=optics.abs_prop)

        return elements


def men_database(database: Path = men_database_path):

    df = read_csv(database)
    # df = df[df['mirror_shape'] != 'Flat']

    return df


def get_men_optics(index_number: int, secondary_contour_points: int = 120):

    men_df = men_database()
    optic_data = men_df[men_df['no'] == index_number]

    # primary field main data #####################################
    number_of_primary_mirrors = 25  # Table 1, p. 594.
    mirror_width = optic_data['W'].item() * 1e3  # in mm
    mirror_center_distance = optic_data['D'].item() * 1e3  # in mm

    if optic_data['mirror_shape'].item() == 'Flat':
        mirror_radius = 0.
    elif optic_data['mirror_shape'].item() == 'Cylindrical':
        mirror_radius = optic_data['MSF_i'].item()
    elif optic_data['mirror_shape'].item() == 'Parabolic':
        mirror_radius = 2 * optic_data['MSF_i'].item()
    else:
        raise ValueError('Undefined shape for the primary mirrors! Please check reference values in the datatable.')

    mirror_radius = mirror_radius * 1e3  # in mm
    ################################################################

    # receiver main data #########################################################################
    receiver_height = optic_data['H_R'].item() * 1e3  # in mm

    absorber_radius = optic_data['r_a'].item() * 1e3  # in mm
    cover_outer_radius = absorber_radius + (0.0225 * 1e3)  # in mm -- Section 3.3, p. 597.
    cover_inner_radius = cover_outer_radius - 3  # in mm -- Table 1, p. 594.

    secondary_optic_gap = optic_data['r'].item() * 1e3  # in mm

    theta_a = optic_data['theta_a'].item()  # in degrees
    theta_max = optic_data['theta_max'].item()  # in degrees
    ###############################################################################################

    tracking_type = optic_data['tracking'].item()

    men_optic = lfc_optic(nbr_mirrors=number_of_primary_mirrors,
                          name=f'men_optic_{index_number}',
                          mirror_width=mirror_width, center_distance=mirror_center_distance,
                          mirror_radius=mirror_radius,
                          receiver_height=receiver_height,
                          absorber_radius=absorber_radius,
                          outer_cover_radius=cover_outer_radius, inner_cover_radius=cover_inner_radius)
    men_optic.tracking_type = tracking_type
    men_optic.index = index_number
    men_optic.add_oommen_cpc(theta_a=theta_a, theta_max=theta_max, gap_radius=secondary_optic_gap,
                             nbr_contour_points=secondary_contour_points)

    # CPC and CEC optics with the same gap as reported by Men ####################
    men_cpc = deepcopy(men_optic)
    men_cpc.name = f'men_{index_number}_cpc'
    men_cpc.add_cpc_secondary(gap_radius=secondary_optic_gap,
                              nbr_contour_points=secondary_contour_points)

    men_cec = deepcopy(men_optic)
    men_cec.name = f'men_{index_number}_cec'
    men_cec.add_cec_secondary(gap_radius=secondary_optic_gap,
                              nbr_contour_points=secondary_contour_points)
    ###############################################################################

    # CPC and CEC with the minimum possible gap #########################################
    # i.e., gap radius equal to the outer cover radius
    # plus the minimum value considered in the interval for $\Delta R_{OD}$ = [0.01, 0.1] (unit in meters)
    # see Table 2, p. 596.
    men_cpc_red_gap = deepcopy(men_optic)
    men_cpc_red_gap.name = f'men_{index_number}_cpc_red_gap'
    men_cpc_red_gap.add_cpc_secondary(gap_radius=cover_outer_radius + 0.01 * 1e3,
                                      nbr_contour_points=secondary_contour_points)

    men_cec_red_gap = deepcopy(men_optic)
    men_cec_red_gap.name = f'men_{index_number}_cec_red_gap'
    men_cec_red_gap.add_cec_secondary(gap_radius=cover_outer_radius + 0.01 * 1e3,
                                      nbr_contour_points=secondary_contour_points)

    ######################################################################################

    return men_optic, men_cpc, men_cpc_red_gap, men_cec, men_cec_red_gap


def simulate_flux(file_name: str, file_path: Path,
                  theta_t: float, dni: float,
                  lfc: lfc_optic,
                  source: RadialSource,
                  optics: OpticalSettings,
                  trace: Trace,
                  rec_dx=0.0, rec_dy=0.0,
                  soltrace_version=2012,
                  force_sim=False) -> Tuple[Path, Path]:

    # Ensuring that argument 'file_path' exist #######
    files_path = Path(file_path)
    files_path.mkdir(parents=True, exist_ok=True)
    ##################################################

    # The Sun #########################################################
    # the sun direction as a [x,y,z] vector
    sun_dir = sun_direction(theta_t=theta_t, theta_l=0.)
    # the object to represent Soltrace Sun box
    sun_box = source.to_soltrace(sun_dir=sun_dir)
    ###################################################################

    # The Optics
    optics_box = Optics(properties=optics.properties)

    # The Geometry
    elements = lfc.to_soltrace(sun_dir=sun_dir, optics=optics, rec_dx=rec_dx, rec_dy=rec_dy)
    stage = Stage(name='linear_fresnel', elements=elements)
    geometry_box = Geometry(stages=[stage])

    # The Element intersections stats
    absorber_stats = ElementStats(stats_name=f'absorber_flux_{file_name}_{theta_t}',
                                  stage_index=0,
                                  element_index=len(elements) - 1,
                                  dni=dni,
                                  x_bins=150,
                                  y_bins=150,
                                  final_rays=True)

    script_path = soltrace_script(sun=sun_box,
                                  optics=optics_box,
                                  geometry=geometry_box,
                                  trace=trace,
                                  stats=absorber_stats,
                                  file_name=f'{file_name}_{theta_t}',
                                  file_path=files_path,
                                  version=soltrace_version)

    # Run SolTrace from the windows prompt cmd ####################################################
    # Only the 2012.7.9 version of SolTrace can run a lk script from the prompt
    # If the absorber_stats file already exist, it should not run SolTrace
    # Finally, if forced, it will always run the simulation
    if (soltrace_version == 2012 and not absorber_stats.file_full_path.is_file()) or force_sim:
        run_soltrace(lk_file_full_path=script_path)
    ################################################################################################

    return script_path, absorber_stats.file_full_path


def optical_analysis(file_name: str, file_path: Path,
                     dni: float,
                     lfc: lfc_optic,
                     source: RadialSource,
                     optics: OpticalSettings,
                     trace: Trace,
                     rec_dx=0., rec_dy=0.,
                     soltrace_version=2012,
                     force_sim=False):

    files_path = Path(file_path, 'flux', lfc.name)
    files_path.mkdir(parents=True, exist_ok=True)

    angle_range = arange(start=0., stop=90., step=5)

    scripts = [0] * angle_range.shape[0]
    stats = [0] * angle_range.shape[0]

    now = datetime.now()
    print(f'It is now {now}, and flux analysis simulations for {lfc.name} has began.')

    for i, theta_t in enumerate(tqdm(angle_range)):

        scripts[i], stats[i] = simulate_flux(file_name=file_name, file_path=files_path,
                                             theta_t=theta_t, dni=dni,
                                             lfc=lfc, rec_dx=rec_dx, rec_dy=rec_dy,
                                             source=source, optics=optics, trace=trace,
                                             soltrace_version=soltrace_version, force_sim=force_sim)

    time.sleep(1.)
    return scripts, stats


def get_optical_analysis_data(file_name: str, file_path: Path,
                              dni: float,
                              lfc: lfc_optic,
                              source: RadialSource,
                              optics: OpticalSettings,
                              trace: Trace,
                              rec_dx=0., rec_dy=0.,
                              soltrace_version=2012):

    files_path = Path(file_path, 'flux', lfc.name)
    files_path.mkdir(parents=True, exist_ok=True)

    angle_range = arange(start=0., stop=90., step=5)

    mirror_area = (lfc.primary_field.widths.sum() * lfc.concentrator_length) / 1e6  # im m2

    efficiency_data = zeros(shape=(angle_range.shape[0], 2))
    uniformity_data = zeros(shape=(angle_range.shape[0], 2))

    for i, theta in enumerate(angle_range):
        absorber_flux_file_path = simulate_flux(file_name=file_name, file_path=files_path,
                                                theta_t=theta, dni=dni,
                                                lfc=lfc, rec_dx=rec_dx, rec_dy=rec_dy,
                                                source=source, optics=optics, trace=trace,
                                                soltrace_version=soltrace_version)[1]

        if not absorber_flux_file_path.is_file():
            raise ValueError('Absorber flux stats does not exist. Please, run simulation before read stats file!')

        absorber_flux_data = ElementFlux(absorber_flux_file_path)

        efficiency_data[i] = theta, absorber_flux_data.total_flux * 1e3 / (mirror_area * dni)
        uniformity_data[i] = theta, absorber_flux_data.flux_x_uniformity

    return efficiency_data, uniformity_data


def simulate_acceptance(file_name: str, file_path: Path,
                        theta_t: float, dni: float,
                        lfc: lfc_optic,
                        source: RadialSource,
                        optics: OpticalSettings,
                        trace: Trace,
                        dt=0.1, lower_value=0.8,
                        soltrace_version=2012,
                        force_sim=False):

    # Ensuring that argument 'file_path' exist #######
    files_path = Path(file_path, f'{theta_t}')
    files_path.mkdir(parents=True, exist_ok=True)
    ##################################################

    # The Optics #######################################
    optics_box = Optics(properties=optics.properties)
    #####################################################

    # The Geometry #######################################################
    sun_on_axis_dir = sun_direction(theta_t=theta_t, theta_l=0.)

    elements = lfc.to_soltrace(sun_dir=sun_on_axis_dir, optics=optics)
    stage = Stage(name='linear_fresnel', elements=elements)
    geometry_box = Geometry(stages=[stage])
    ######################################################################

    # the sun direction as a [x,y,z] vector
    # the object to represent Soltrace Sun box
    ###################################################################

    on_axis_stats_file = simulate_flux(file_name=file_name, file_path=files_path,
                                       lfc=lfc,
                                       theta_t=theta_t, dni=dni,
                                       source=source, optics=optics, trace=trace)[1]
    on_axis_flux = ElementFlux(stats_file=on_axis_stats_file).total_flux

    angles = [0]
    norm_flux = [1]

    # loop for the positive off-axis incidence variation
    k = 1
    while norm_flux[-1] >= lower_value:

        # off_axis incidence
        off_axis_dev = round(k * dt, 5)
        inc_angle = theta_t + off_axis_dev
        sun_dir = sun_direction(theta_t=inc_angle, theta_l=0.)
        sun_box = source.to_soltrace(sun_dir=sun_dir)
        # the tracking angle of the primary field does not change

        # The Element stats
        absorber_stats = ElementStats(stats_name=f'absorber_flux_{file_name}_{off_axis_dev}',
                                      stage_index=0,
                                      element_index=len(elements) - 1,
                                      dni=dni,
                                      x_bins=150,
                                      y_bins=150,
                                      final_rays=True)

        script_path = soltrace_script(sun=sun_box,
                                      optics=optics_box,
                                      geometry=geometry_box,
                                      trace=trace,
                                      stats=absorber_stats,
                                      file_name=f'{file_name}_{off_axis_dev}',
                                      file_path=files_path,
                                      version=soltrace_version)

        if (soltrace_version == 2012 and not absorber_stats.file_full_path.is_file()) or force_sim:
            run_soltrace(lk_file_full_path=script_path)

        # getting the total flux at the absorber for the off-axis incidence
        absorber_flux = ElementFlux(stats_file=absorber_stats.file_full_path).total_flux

        # updating the values of incidence angle and normalized flux at the absorber
        angles.append(k * dt)
        norm_flux.append(absorber_flux / on_axis_flux)
        k += 1

    ################################################################################################

    # loop for the negative off-axis incidence variation
    k = 1
    while norm_flux[0] >= lower_value:

        # off_axis incidence
        off_axis_dev = -round(k * dt, 5)
        inc_angle = theta_t + off_axis_dev
        sun_dir = sun_direction(theta_t=inc_angle, theta_l=0.)
        sun_box = source.to_soltrace(sun_dir=sun_dir)
        # the tracking angle of the primary field does not change

        # The Element stats
        absorber_stats = ElementStats(stats_name=f'absorber_flux_{file_name}_{off_axis_dev}',
                                      stage_index=0,
                                      element_index=len(elements) - 1,
                                      dni=dni,
                                      x_bins=150,
                                      y_bins=150,
                                      final_rays=True)

        script_path = soltrace_script(sun=sun_box,
                                      optics=optics_box,
                                      geometry=geometry_box,
                                      trace=trace,
                                      stats=absorber_stats,
                                      file_name=f'{file_name}_{off_axis_dev}',
                                      file_path=files_path,
                                      version=soltrace_version)

        if (soltrace_version == 2012 and not absorber_stats.file_full_path.is_file()) or force_sim:
            run_soltrace(lk_file_full_path=script_path)

        # getting the total flux at the absorber for the off-axis incidence
        absorber_flux = ElementFlux(stats_file=absorber_stats.file_full_path).total_flux

        # updating the values of incidence angle and normalized flux at the absorber
        angles.insert(0, - k * dt)
        norm_flux.insert(0, absorber_flux / on_axis_flux)

        # updating loop control variable
        k += 1

    # Returning simulation results
    transmission_data = zeros(shape=(len(angles), 2))
    transmission_data.T[:] = angles, norm_flux

    return transmission_data


def acceptance_analysis(file_name: str, file_path: Path,
                        dni: float,
                        lfc: lfc_optic,
                        source: RadialSource,
                        optics: OpticalSettings,
                        trace: Trace,
                        dt=0.1, lower_value=0.8,
                        soltrace_version=2012,
                        force_sim=False) -> array:

    simulation_files_path = Path(file_path, 'acceptance', lfc.name)

    angle_range = arange(start=0., stop=90., step=5)

    acceptance_half_angles = zeros(shape=(angle_range.shape[0], 2))

    now = datetime.now()
    print(f'It is now {now}, and acceptance analysis simulations for {lfc.name} has began.')

    for i, theta in enumerate(tqdm(angle_range)):

        acceptance_data = simulate_acceptance(lfc=lfc,
                                              theta_t=theta, dni=dni,
                                              file_name=file_name, file_path=simulation_files_path,
                                              source=source, optics=optics, trace=trace,
                                              dt=dt, lower_value=lower_value,
                                              soltrace_version=soltrace_version, force_sim=force_sim)

        theta_a = acceptance_angle(acceptance_data=acceptance_data, ref_value=0.9)

        acceptance_half_angles[i] = theta, theta_a

    time.sleep(1.)

    return acceptance_half_angles


def get_simulations_data(file_name: str, file_path: Path,
                         dni: float,
                         lfc: lfc_optic,
                         source: RadialSource,
                         optics: OpticalSettings,
                         trace: Trace,
                         angle_step: float, lower_value: float,
                         soltrace_version=2012):

    efficiency_data, uniformity_data = get_optical_analysis_data(lfc=lfc, dni=dni,
                                                                 file_name=file_name, file_path=file_path,
                                                                 source=source, optics=optics,
                                                                 trace=trace, soltrace_version=soltrace_version)

    acceptance_data = acceptance_analysis(file_name=file_name, file_path=file_path,
                                          dni=dni, lfc=lfc, source=source, optics=optics,
                                          trace=trace,
                                          dt=angle_step, lower_value=lower_value)

    return efficiency_data, uniformity_data, acceptance_data


def get_averaged_values(site_data: SiteData,
                        file_name: str, file_path: Path,
                        dni: float,
                        lfc: lfc_optic,
                        source: RadialSource,
                        optics: OpticalSettings,
                        trace: Trace,
                        angle_step: float, lower_value: float,
                        soltrace_version=2012):

    # Getting performance data
    efficiency_data, uniformity_data, acceptance_data = get_simulations_data(file_name=file_name, file_path=file_path,
                                                                             dni=dni, lfc=lfc, source=source,
                                                                             optics=optics, trace=trace,
                                                                             angle_step=angle_step,
                                                                             lower_value=lower_value,
                                                                             soltrace_version=soltrace_version)

    # Creating the interpolation function for each performance metric
    eff_function = interp1d(*efficiency_data.T, kind='linear')
    uni_function = interp1d(*uniformity_data.T, kind='linear')
    acc_function = interp1d(*acceptance_data.T, kind='linear')

    # the transversal and longitudinal incidence angles for the location
    transversal_angles, longitudinal_angles = site_data.linear_angles(NS=True, solar_longitudinal=False)

    # adding the linear angles to the location tmy
    # concentrator is symmetric
    tmy = site_data.tmy_data
    tmy['theta_t'] = absolute(transversal_angles)
    tmy['theta_l'] = absolute(longitudinal_angles)

    # ns_calculations
    ns_tmy = tmy[tmy['theta_t'] <= 85.]
    angles = ns_tmy['theta_t']
    dni = ns_tmy['dni'].array

    ns_eff = eff_function(angles).dot(dni) / dni.sum()
    ns_uni = uni_function(angles).dot(dni) / dni.sum()
    ns_acc = acc_function(angles).dot(dni) / dni.sum()

    # ew_calculations
    ew_tmy = tmy[tmy['theta_l'] <= 85.]
    angles = ew_tmy['theta_l']
    dni = ew_tmy['dni'].array

    ew_eff = eff_function(angles).dot(dni) / dni.sum()
    ew_uni = uni_function(angles).dot(dni) / dni.sum()
    ew_acc = acc_function(angles).dot(dni) / dni.sum()

    # output dictionaries

    ns_dic = {'avg_eff': ns_eff, 'avg_uni': ns_uni, 'avg_acc': ns_acc}
    ew_dic = {'avg_eff': ew_eff, 'avg_uni': ew_uni, 'avg_acc': ew_acc}

    return ns_dic, ew_dic
