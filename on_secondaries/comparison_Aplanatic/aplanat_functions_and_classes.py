import time
from datetime import datetime
from pathlib import Path
from typing import Tuple

from numpy import arange, absolute
from scipy.interpolate import interp1d
from tqdm import tqdm

from solaren.niopy.geometric_transforms import dst, Sy, isl, tgs2tube
from solaren.niopy.plane_curves import concatenate_curves, PlaneCurve
from solaren.scopy import OpticalProperty
from solaren.scopy.aplanats import *
from solaren.scopy.linear_fresnel import Absorber, rabl_curvature, PrimaryField, acceptance_angle, PrimaryMirror
from solaren.scopy.linear_fresnel.secondaries import half_acceptance2tube, cpc4tube, cec4tube
from solaren.scopy.sunlight import RadialSource, sun_direction, SiteData
from solaren.soltracepy import Optics, Stage, Geometry, ElementStats, soltrace_script, run_soltrace, ElementFlux, Trace


class OpticalSettings:

    def __init__(self,
                 primaries_property: OpticalProperty.reflector,
                 secondary_property: OpticalProperty.secondary,
                 absorber_property: OpticalProperty.absorber_tube):

        self.primaries_property = primaries_property
        self.secondary_property = secondary_property
        self.absorber_property = absorber_property

        self.primaries_prop = primaries_property.to_soltrace()
        self.secondaries_prop = secondary_property.to_soltrace()
        self.abs_prop = self.absorber_property.to_soltrace()

        self.properties = [self.primaries_prop,
                           self.secondaries_prop,
                           self.abs_prop]


class AplanatOptic:

    def __init__(self,
                 s: float,
                 k: float,
                 NA: float,
                 primary_width: float, contour_points=120):

        # Attributes to hold initializer arguments ####
        self.s = s
        self.k = k
        self.na = abs(NA)
        self.primary_width = primary_width
        #########################################################

        # Auxiliary attributes ##################################
        self.phi_max = arcsin(self.na)
        self.phi_values = linspace(start=0., stop=self.phi_max, num=contour_points)

        self.xp_max = primary_width / 2.
        self.f = self.xp_max / sin(self.phi_max)
        #########################################################

        # the contour of aplanat primary and secondary optics
        self.primary_contour, self.secondary_contour = self.dual4tube(contour_points=contour_points)

    def aplanatic_functions(self):
        return dual_aplanat_functions_for_tube(s=self.s, k=self.k)

    def dual4tube(self, contour_points: int = 150):

        half_contour_points = contour_points // 2
        right_primary, left_secondary = dual4tube(s=self.s,
                                                  k=self.k,
                                                  NA=self.na,
                                                  nbr_pts=half_contour_points)

        right_secondary = array([Sy(pt) for pt in left_secondary])
        left_primary = array([Sy(pt) for pt in right_primary])

        primary_contour = concatenate_curves(base_curve=left_primary[::-1], next_curve=right_primary)
        secondary_contour = concatenate_curves(base_curve=left_secondary[::-1], next_curve=right_secondary)

        return primary_contour * self.f, secondary_contour * self.f

    def segment_primary(self, nbr_segments: int):

        # Primary field values in the x-axis ########
        x_min = self.primary_contour.T[0].min()
        x_max = self.primary_contour.T[0].max()
        #############################################

        # calculating the size of the segment in the x-axis
        dx = (x_max - x_min) / nbr_segments

        # x_values of the center of the segments
        seg_x_values = array([(x_min + 0.5*dx) + i * dx for i in range(nbr_segments)])

        # the value of the parameter 'phi' for those centers
        seg_phi_values = arcsin(seg_x_values / self.f)[::-1]

        # call for the base aplanat unit equations
        primary_function, secondary_function = self.aplanatic_functions()

        # Calculates the center of the segments in the primary optic contour -- and scale it
        segments_centers = array([primary_function(phi) for phi in seg_phi_values]) * self.f

        # calculating the corresponding points in the secondary optic contour -- and scale it
        tracking_points = array([secondary_function(phi) for phi in seg_phi_values]) * self.f

        return segments_centers, tracking_points

    def linear_fresnel(self, number_primaries: int):

        # The height difference between the focus and the lowest point in the aplanat primary optic
        h = abs(self.k - self.s) * self.f

        # segment the primary optic in a finite number of segments
        segments_centers, tracking_points = self.segment_primary(nbr_segments=number_primaries)

        # The lowest point in the aplanat primary optic
        p = array([0, -h])
        # a vector to define a horizontal line in which the centers of the primary mirrors will be defined
        v = array([1, 0])

        # the center points of the primaries are defined as the interception between two straight lines
        # one is the line defined by the horizontal line
        # the other is the one which connects the center of the segment with the corresponding point in the secondary
        horizontal_centers = array([isl(p=p, v=v, q=c, u=c - a)
                                    for c, a in zip(segments_centers, tracking_points)])

        # Creating the primary field
        primaries_width = 2 * (self.xp_max - segments_centers.T[0].max())
        primaries_radii = [rabl_curvature(center=c, aim=a, theta_d=0)
                           for c, a in zip(horizontal_centers, tracking_points)]

        heliostats = [PrimaryMirror(width=primaries_width,
                                    center=c, radius=r) for c, r in zip(horizontal_centers, primaries_radii)]

        primary_field = PrimaryField(heliostats=heliostats)

        return tracking_points, primary_field


class AplanatLfc:

    def __init__(self,
                 name: str,
                 s: float,
                 k: float,
                 NA: float,
                 aplanat_primary_width: float,
                 absorber_radius: float,
                 nbr_primary_mirrors: int,
                 contour_points: int,
                 concentrator_length: float = 6000):

        # auxiliary attributes
        self.name = name
        self.concentrator_length = concentrator_length

        # The base AplanatOptic object ############################################
        self.base_aplanat = AplanatOptic(s=s, k=k, NA=NA,
                                         primary_width=aplanat_primary_width,
                                         contour_points=contour_points)
        ###########################################################################

        # An object for the absorber tube ##########################################
        self.absorber = Absorber.tube(radius=absorber_radius,
                                      center=array([0, 0]))
        ############################################################################

        # The primary field and the corresponding tracking points ######################################

        # tracking points in the secondary contour and the primary field object
        aplanat_lfc_data = self.base_aplanat.linear_fresnel(number_primaries=nbr_primary_mirrors)
        self.tracking_points, self.primary_field = aplanat_lfc_data

        # the primary field total width and edge-points
        self.primary_field_width = dst(p=self.primary_field.centers[0], q=self.primary_field.centers[-1])
        self.primary_field_width += self.primary_field.widths.mean()

        self.f1 = array([-0.5 * self.primary_field_width, self.primary_field.centers[0][-1]])
        self.f2 = array([+0.5 * self.primary_field_width, self.primary_field.centers[0][-1]])
        ################################################################################################

        # The secondary optic ####################################################################################
        self.secondary_type = 'aplanat'
        self.secondary_center = self.absorber.center - array([0, self.base_aplanat.k * self.base_aplanat.f])
        self.secondary = PlaneCurve(curve_pts=self.base_aplanat.secondary_contour,
                                    curve_center=self.secondary_center)
        ##########################################################################################################

    def tube_edges(self):

        # Right edge of the receiver  ##############################################
        # it considers the left edge of the primary field
        ta, tb = tgs2tube(point=self.f1,
                          tube_radius=self.absorber.radius,
                          tube_center=self.absorber.center)

        if ta[0] < tb[0]:
            t2 = tb
        else:
            t2 = ta
        ############################################################################

        # Left edge of the receiver #################################################
        # it considers the right edge of the primary field

        ta, tb = tgs2tube(point=self.f2,
                          tube_radius=self.absorber.radius,
                          tube_center=self.absorber.center)

        if ta[0] < tb[0]:
            t1 = ta
        else:
            t1 = tb
        ############################################################################

        return t1, t2

    def edge_rays_acceptance(self) -> float:

        half_acceptance = half_acceptance2tube(primary_field=self.primary_field.primaries,
                                               tube_center=self.absorber.center,
                                               tube_radius=self.absorber.radius)

        return half_acceptance

    def add_cpc_secondary(self,
                          name: str = None,
                          gap_radius: float = None, nbr_contour_points: int = 120):

        # Secondary optic design ##########################################################################
        # cpc secondary is design with four sections: two involutes and two macro focal parabolas
        points_per_section = int(round(nbr_contour_points / 4))

        # checking if the gap is forced or should follow the highest point in the aplanat secondary optic
        g = dst(self.secondary_center, self.absorber.center) if gap_radius is None else gap_radius

        # design the CPC secondary optic
        secondary_curve = cpc4tube(primary_field=self.primary_field.primaries,
                                   tube_center=self.absorber.center,
                                   tube_radius=self.absorber.radius,
                                   gap_radius=g,
                                   points_per_section=points_per_section)
        ######################################################################################################

        # Updating the attributes for the addition of a new secondary optic ##################################
        self.name = 'aplanat_cpc' if name is None else name

        self.secondary_type = 'cpc'
        self.tracking_points = [secondary_curve.aperture_center.round(10)] * len(self.primary_field.centers)
        self.secondary = secondary_curve.as_plane_curve()
        ######################################################################################################

        # Updating the primary field ###########################################################################

        primaries_radii = [rabl_curvature(center=c, aim=a, theta_d=0)
                           for c, a in zip(self.primary_field.centers, self.tracking_points)]

        heliostats = [PrimaryMirror(width=self.primary_field.widths.mean(),
                                    center=c, radius=r) for c, r in zip(self.primary_field.centers, primaries_radii)]

        self.primary_field = PrimaryField(heliostats=heliostats)
        #######################################################################################################

        return self.secondary

    def add_cec_secondary(self,
                          name: str = None,
                          gap_radius: float = None, nbr_contour_points: int = 120):

        # Secondary optic design #############################################################################
        # cec secondary is design with four sections: two involutes and two macro focal ellipses
        points_per_section = int(round(nbr_contour_points / 4))

        # checking if the gap is forced or should follow the highest point in the aplanat secondary optic
        g = dst(self.secondary_center, self.absorber.center) if gap_radius is None else gap_radius

        # design the CEC secondary optic
        secondary_curve = cec4tube(primary_field=self.primary_field.primaries,
                                   tube_center=self.absorber.center,
                                   tube_radius=self.absorber.radius,
                                   gap_radius=g,
                                   points_per_section=points_per_section)
        ######################################################################################################

        # Updating the attributes for the addition of a new secondary optic ##################################
        self.name = 'aplanat_cec' if name is None else name

        self.secondary_type = 'cec'
        self.tracking_points = [secondary_curve.aperture_center.round(10)] * len(self.primary_field.centers)
        self.secondary = secondary_curve.as_plane_curve()
        ######################################################################################################

        # Updating the primary field ###########################################################################

        primaries_radii = [rabl_curvature(center=c, aim=a, theta_d=0)
                           for c, a in zip(self.primary_field.centers, self.tracking_points)]

        heliostats = [PrimaryMirror(width=self.primary_field.widths.mean(),
                                    center=c, radius=r) for c, r in zip(self.primary_field.centers, primaries_radii)]

        self.primary_field = PrimaryField(heliostats=heliostats)
        #######################################################################################################

        return self.secondary

    def to_soltrace(self, sun_dir: array, optics: OpticalSettings,
                    rec_dx=0., rec_dy=0.) -> list:

        # The primary mirrors as SolTrace elements
        elements = self.primary_field.to_soltrace(rec_aim=self.tracking_points,
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
                                                          optic=optics.abs_prop,
                                                          name='absorber_tube')
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

            displaced_absorber = Absorber.tube(radius=self.absorber.radius,
                                               center=self.absorber.center + array([rec_dx, rec_dy]))

            elements += displaced_absorber.as_soltrace_element(length=self.concentrator_length,
                                                               optic=optics.abs_prop,
                                                               name='absorber_tube')

        return elements


def simulate_flux(file_name: str, file_path: Path,
                  theta_t: float, dni: float,
                  lfc: AplanatLfc,
                  source: RadialSource,
                  optics: OpticalSettings,
                  trace: Trace,
                  rec_dx=0., rec_dy=0.,
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
                     lfc: AplanatLfc,
                     source: RadialSource,
                     optics: OpticalSettings,
                     trace: Trace,
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
                                             lfc=lfc,
                                             source=source, optics=optics, trace=trace,
                                             soltrace_version=soltrace_version, force_sim=force_sim)

    time.sleep(1.)
    return scripts, stats


def get_optical_analysis_data(file_name: str, file_path: Path,
                              dni: float,
                              lfc: AplanatLfc,
                              source: RadialSource,
                              optics: OpticalSettings,
                              trace: Trace,
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
                                                lfc=lfc,
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
                        lfc: AplanatLfc,
                        source: RadialSource,
                        optics: OpticalSettings,
                        trace: Trace,
                        dt=0.1, lower_value=0.8,
                        soltrace_version=2012,
                        force_sim=False):

    # Ensuring that argument 'file_path' exist #######
    files_path = Path(file_path)
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
        inc_angle = theta_t + k * dt
        sun_dir = sun_direction(theta_t=inc_angle, theta_l=0.)
        sun_box = source.to_soltrace(sun_dir=sun_dir)
        # the tracking angle of the primary field does not change

        # The Element stats
        absorber_stats = ElementStats(stats_name=f'absorber_flux_{file_name}_{theta_t}_{k * dt}',
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
                                      file_name=f'{file_name}_{theta_t}_{k * dt}',
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
        inc_angle = theta_t - k * dt
        sun_dir = sun_direction(theta_t=inc_angle, theta_l=0.)
        sun_box = source.to_soltrace(sun_dir=sun_dir)
        # the tracking angle of the primary field does not change

        # The Element stats
        absorber_stats = ElementStats(stats_name=f'absorber_flux_{file_name}_{theta_t}_{-k * dt}',
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
                                      file_name=f'{file_name}_{theta_t}_{-k * dt}',
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
                        lfc: AplanatLfc,
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
                         lfc: AplanatLfc,
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
                        lfc: AplanatLfc,
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
    ns_tmy = tmy[tmy['theta_t'] <= 85]
    angles = ns_tmy['theta_t']
    dni = ns_tmy['dni']

    ns_eff = eff_function(angles).dot(dni) / dni.sum()
    ns_uni = uni_function(angles).dot(dni) / dni.sum()
    ns_acc = acc_function(angles).dot(dni) / dni.sum()

    # ew_calculations
    ew_tmy = tmy[tmy['theta_l'] <= 85]
    angles = ew_tmy['theta_l']
    dni = ew_tmy['dni']

    ew_eff = eff_function(angles).dot(dni) / dni.sum()
    ew_uni = uni_function(angles).dot(dni) / dni.sum()
    ew_acc = acc_function(angles).dot(dni) / dni.sum()

    # output dictionaries

    ns_dic = {'avg_eff': ns_eff, 'avg_uni': ns_uni, 'avg_acc': ns_acc}
    ew_dic = {'avg_eff': ew_eff, 'avg_uni': ew_uni, 'avg_acc': ew_acc}

    return ns_dic, ew_dic
