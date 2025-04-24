import json
from datetime import datetime
from pathlib import Path

from numpy import array, ones, arange, zeros, ndarray, linspace, deg2rad
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import fsolve

from tqdm import tqdm

from solaren.niopy.geometric_transforms import dst
from solaren.scopy.linear_fresnel import (Absorber, uniform_centers, LFR, rabl_curvature, PrimaryMirror,
                                          PrimaryField, annual_eta, boito_curvature)

from solaren.scopy.sunlight import RadialSource, SiteData
from utils import dic2json

########################################################################################################################
########################################################################################################################


# Classes and functions for this program ###############################################################################


class lfr_geometry:
    """
    This class represents the main data of a geometrical configuration regarding the curvature design.
    In this sense, it also encompasses the primary field geometric data and the flat, horizontal receiver.

    It considers a primary field with uniform configurations of width and shift (distance between the center points of
    two neighbor mirrors).

    """

    def __init__(self, name: str, rec_width: float, rec_center: array,
                 mirror_width: float, nbr_mirrors: int,
                 total_width: float = None, center_distance: float = None):
        # The name of the linear fresnel geometry ####
        self.name = name
        ##############################################

        # receiver data and object #############################
        self.rec_width = abs(rec_width)
        self.rec_center = rec_center

        self.receiver = Absorber.flat(width=rec_width,
                                      center=rec_center)
        self.rec_aim = self.receiver.sm
        ########################################################

        # primary field main data ################################################################################
        self.mirror_width = abs(mirror_width)
        self.nbr_mirrors = abs(nbr_mirrors)
        self.widths = ones(self.nbr_mirrors) * self.mirror_width

        # calculating the center points of the primary field mirrors ##########
        self.centers = uniform_centers(mirror_width=self.mirror_width,
                                       nbr_mirrors=self.nbr_mirrors,
                                       total_width=total_width,
                                       center_distance=center_distance)

        # calculating the correspondent shift and total width of the primary field
        self.center_distance = dst(p=self.centers[0], q=self.centers[1])
        self.total_width = dst(p=self.centers[0], q=self.centers[-1]) + self.mirror_width

        # filling factor of the primary field - the ratio between reflective and occupied aperture
        self.filling_factor = self.widths.sum() / self.total_width
        ##########################################################################################################

        # the edge-points of the primary field ############
        self.f1 = array([-0.5 * self.total_width, 0])
        self.f2 = array([+0.5 * self.total_width, 0])
        ###################################################

    # A method to export geometry as a dictionary

    def export2dic(self):
        dic = {'name': self.name,
               'mirror_width': self.mirror_width,
               'nbr_mirrors': self.nbr_mirrors,
               'total_width': self.total_width,
               'rec_width': self.rec_width,
               'rec_height': self.rec_center[-1],
               'units': 'millimeters'}

        return dic

    def export2json(self, file_path):
        file_full_path = Path(file_path, f'{self.name}_geometry.json')

        dic = self.export2dic()
        dic2json(d=dic, file_path=file_path, file_name=f'{self.name}_geometry')

        return file_full_path


def design_position_curvature(geom: lfr_geometry, curvature_design) -> LFR:
    """
    This function calculates the curvature radius of each mirror in a primary field for a particular curvature design,
    and then it returns the correspondent LFR object (primary field and flat receiver).

    :param geom: a linear fresnel geometry, a lfr_geometry object
    :param curvature_design: the curvature design.

    :return: an LFR object
    """

    # list of heliostat to create the primary field
    heliostats = []

    # variables to design the Heliostat objects
    centers = geom.centers
    widths = geom.widths
    sm = geom.receiver.sm

    # Creating the Heliostat objects
    for w, hc in zip(widths, centers):

        if isinstance(curvature_design, (float, int)):
            rr = rabl_curvature(center=hc,
                                aim=sm,
                                theta_d=curvature_design)
        elif curvature_design == 'SR':
            rr = 2 * dst(hc, sm)
        elif curvature_design == 'flat':
            rr = 0.

        else:
            raise ValueError("Design position must be an number, 'SR' or 'flat'.")

        # append the current heliostat to the list
        hel = PrimaryMirror(center=hc,
                            width=w,
                            radius=rr)
        heliostats.append(hel)

    # PrimaryField and LFR instances
    primary_field = PrimaryField(heliostats=heliostats)
    lfr_concentrator = LFR(primary_field=primary_field, flat_absorber=geom.receiver)

    return lfr_concentrator


def design_positions_analysis(geom: lfr_geometry, site: SiteData, cum_eff: array):
    """
    This function calculates the annual efficiency for all possible curvature designs regarding a design position in
    the range [-85º, +85º] for both North-South and East-West oriented fields.

    :param geom: a linear fresnel geometry.
    :param site: a SiteData object with weather data of a particular location.
    :param cum_eff: a cumulative density function of the considered effective source.

    :return: annual efficiency data for NS and EW mountings.
    """

    # range of designs positions to be analyzed
    theta_i, theta_f, dt = -85., 85., 2.5
    theta_d_range = arange(start=theta_i, stop=theta_f + dt, step=dt)

    # computation of the factorized efficiencies of these curvature designs
    fac_etas = [design_position_curvature(geom=geom,
                                          curvature_design=td).optical_analysis(cum_eff=cum_eff,
                                                                                rec_aim=geom.rec_aim,
                                                                                end_losses=False,
                                                                                factorized=True)
                for td in tqdm(theta_d_range)]

    # calculation of the NS annual efficiencies
    ns_data = zeros(shape=(theta_d_range.shape[0], 2))
    ns_etas = [annual_eta(transversal_data=t[0], longitudinal_data=t[1], site=site, NS=True) for t in fac_etas]
    ns_data.T[:] = theta_d_range, ns_etas

    # calculation of the EW annual efficiencies
    ew_data = zeros(shape=(theta_d_range.shape[0], 2))
    ew_etas = [annual_eta(transversal_data=t[0], longitudinal_data=t[1], site=site, NS=False) for t in fac_etas]
    ew_data.T[:] = theta_d_range, ew_etas

    return ns_data, ew_data


def boito_design(geom: lfr_geometry, lat: float):
    """
    This function calculates the correspondent curvature radius design as proposed by Boito and Grena (2017).

    :param geom: a linear fresnel geometry
    :param lat: local latitude, in radians

    :return: an LFR object
    """

    # Heliostat object list
    heliostats = []

    # calculate curvature radius as proposed by Boito and Grena and design the Heliostat object
    for w, hc in zip(geom.widths, geom.centers):
        r = boito_curvature(center=hc, aim=geom.rec_aim, lat=lat)
        hel = PrimaryMirror(center=hc,
                            width=w,
                            radius=r)
        heliostats.append(hel)

    # LFR object
    primary_field = PrimaryField(heliostats=heliostats)
    lfr = LFR(primary_field=primary_field, flat_absorber=geom.receiver)

    return lfr


def lfr_design(geom: lfr_geometry, radii: array):
    """
    This function designs an LFR considering a particular geometry and a set of curvature radius for the primary field.

    :param geom: a linear fresnel geometry
    :param radii: the list of curvature radius for the primary mirrors

    :return: an LFR object
    """

    widths = geom.widths
    centers = geom.centers

    heliostats = []
    if isinstance(radii, ndarray):
        # Check if the number of primary mirrors matches with the number of curvature radiuses
        assert geom.centers.shape[0] == radii.shape[0], 'The number of mirror does not match with the number of radii!'

        for w, hc, r in zip(widths, centers, radii):
            hel = (PrimaryMirror(width=w,
                                 center=hc,
                                 radius=r))
            heliostats.append(hel)
    elif isinstance(radii, (float, int)):
        for w, hc in zip(widths, centers):
            hel = PrimaryMirror(width=w,
                                center=hc,
                                radius=radii)
            heliostats.append(hel)
    else:
        raise ValueError('Radii argument should be an array, float, or int.')

    primary_field = PrimaryField(heliostats)
    lfr_concentrator = LFR(primary_field=primary_field, flat_absorber=geom.receiver)

    return lfr_concentrator


def heliostat_optimum_radius(i: int, geom: lfr_geometry, site: SiteData, cum_eff: array, init_radii: array = None):
    """
    This functions runs a univariate parametric analysis over the curvature radius of a particular heliostat to
    determine the ones that yield the higher annual optical efficiency for NS and EW mountings.

    :param i: index of the mirror in the primary field
    :param geom: a linear fresnel geometry
    :param site: a SiteData object
    :param cum_eff: a cumulative density function of an effective source
    :param init_radii: an initial set of curvature radius for the primary field

    :return: a tuple of floats with the optimal curvature radius (NS, EW).
    """

    # The lfr curvature radius of all mirrors
    lfr_radii = array([2 * dst(hc, geom.rec_aim) for hc in geom.centers]) if init_radii is None else init_radii

    # the range of curvature radius of the i-th heliostat to be analyzed
    hel_radii = linspace(start=0.75 * lfr_radii[i], stop=3 * lfr_radii[i], num=40)

    # array to hold data of annual efficiency for NS and EW mountings
    ns_etas = zeros(hel_radii.shape[0])
    ew_etas = zeros(hel_radii.shape[0])

    # Univariate parametric analysis
    for j, r in enumerate(hel_radii):
        # the curvature radius of the i-th heliostat being analyzed
        lfr_radii[i] = r

        # the correspondent LFR object
        lfr = lfr_design(geom=geom, radii=lfr_radii)

        # its factorized efficiencies
        fac_eta = lfr.optical_analysis(cum_eff=cum_eff, symmetric=False, factorized=True, end_losses=False)

        # NS and EW annual efficiencies calculations
        ns_etas[j] = annual_eta(transversal_data=fac_eta[0], longitudinal_data=fac_eta[1], site=site, NS=True)
        ew_etas[j] = annual_eta(transversal_data=fac_eta[0], longitudinal_data=fac_eta[1], site=site, NS=False)

    # determining the optical curvature radius for the NS mounting
    ns_etas_function = InterpolatedUnivariateSpline(hel_radii, ns_etas, k=4)
    ns_eta_roots = ns_etas_function.derivative().roots()
    ns_opt_rad = ns_eta_roots[0] if ns_eta_roots.size > 0 else fsolve(ns_etas_function.derivative(), hel_radii[0])[0]

    # determining the optical curvature radius for the EW mounting
    ew_etas_function = InterpolatedUnivariateSpline(hel_radii, ew_etas, k=4)
    ew_eta_roots = ew_etas_function.derivative().roots()
    ew_opt_rad = ew_eta_roots[0] if ew_eta_roots.size > 0 else fsolve(ew_etas_function.derivative(), hel_radii[0])[0]

    return ns_opt_rad, ew_opt_rad


def non_uniform_optimization_routine(geom: lfr_geometry, site: SiteData, cum_eff: array):
    """
    This function computes the non-uniform optimization routine for the non-uniform curvature radius design.
    It runs the function 'heliostat_optimum_radius' for each mirror in the primary field

    :param geom: a linear fresnel geometry
    :param site: a SiteData object
    :param cum_eff: a cumulative density function of an effective source

    :return: a tuple of LFR objects (NS, EW)
    """

    # variables to hold optimum curvature radius data
    ns_opt_radii = zeros(geom.centers.shape[0])
    ew_opt_radii = zeros(geom.centers.shape[0])
    ########################################################

    # run the optimization for each heliostat in the primary field
    for i, hc in enumerate(tqdm(geom.centers)):
        ns_opt_radii[i], ew_opt_radii[i] = heliostat_optimum_radius(i=i, geom=geom, site=site, cum_eff=cum_eff)

    # creating the Heliostat objects with the optimized curvature radius calculated before
    ns_heliostats = []
    ew_heliostats = []
    for w, hc, ns_r, ew_r in zip(geom.widths, geom.centers, ns_opt_radii, ew_opt_radii):
        ns_hel = PrimaryMirror(center=hc, width=w, radius=ns_r)
        ns_heliostats.append(ns_hel)

        ew_hel = PrimaryMirror(center=hc, width=w, radius=ew_r)
        ew_heliostats.append(ew_hel)
    ###############################################################################################

    # creating the LFR objects for the NS and EW cases #########################
    ns_primary_field = PrimaryField(heliostats=ns_heliostats)
    ns_lfr = LFR(primary_field=ns_primary_field, flat_absorber=geom.receiver)

    ew_primary_field = PrimaryField(heliostats=ew_heliostats)
    ew_lfr = LFR(primary_field=ew_primary_field, flat_absorber=geom.receiver)
    ############################################################################

    return ns_lfr, ew_lfr


def non_uniform_design_analysis(geom: lfr_geometry, site: SiteData, source: RadialSource):
    """
    This function computes the non-uniform design analysis, i.e., it computes the annual efficiencies of all possible
    non-uniform curvature designs:
        * Design position range
        * Boito and Grena design
        * Specific reference design
        * Non-uniform optimization routine

    :param geom: a linear fresnel geometry
    :param site: a SiteData object
    :param source: a RadialSource object

    :return: a tuple of dictionaries with (NS, EW) data
    """

    files_path = Path(Path.cwd(), site.name, geom.name, source.name, 'efficiency')
    files_path.mkdir(parents=True, exist_ok=True)

    # the liner cumulative density function of the effective source to perform calculations
    cum_eff = source.linear_cum_eff()

    # NS and EW design positions analysis ############################################
    # Each one of this variables is an array whose shape is (#,2).
    # Each datapoint of these arrays is of the kind [theta_d, eta(theta_d)].

    dp_file = Path(files_path, 'dp_data.json')
    if dp_file.is_file():
        with open(dp_file, 'r') as loaded_file:
            dp_dic = json.load(loaded_file)
            ns_dp = array(dp_dic['ns'])
            ew_dp = array(dp_dic['ew'])
    else:
        ns_dp, ew_dp = design_positions_analysis(geom=geom, site=site, cum_eff=cum_eff)
        dp_dic = {'ns': ns_dp, 'ew': ew_dp}
        dp_dic = dic2json(d=dp_dic, file_path=files_path, file_name='dp_data')
        ns_dp = array(dp_dic['ns'])
        ew_dp = array(dp_dic['ew'])

    # NS optimum design position
    ns_function = InterpolatedUnivariateSpline(*ns_dp.T, k=4)
    ns_opt_dp = ns_function.derivative().roots()

    ns_max_eta = ns_function(ns_opt_dp).max()
    ns_max_dp = ns_opt_dp[ns_function(ns_opt_dp).tolist().index(ns_max_eta)]

    # EW optimum design position

    ew_function = InterpolatedUnivariateSpline(*ew_dp.T, k=4)
    ew_opt_dp = ew_function.derivative().roots()

    ew_max_eta = ew_function(ew_opt_dp).max()
    ew_max_dp = ew_opt_dp[ew_function(ew_opt_dp).tolist().index(ew_max_eta)]

    ##################################################################################

    # Sun reference ###############################

    ns_sun_dp, ew_sun_dp = site.sun_references_positions()

    ns_sun_eta = ns_function(ns_sun_dp)
    ew_sun_eta = ew_function(ew_sun_dp)

    ##################################################################################

    # Boito design #############################################################################################
    boito_fac_etas = boito_design(geom=geom, lat=deg2rad(site.lat)).optical_analysis(cum_eff=cum_eff,
                                                                                     factorized=True,
                                                                                     symmetric=False,
                                                                                     end_losses=False)

    ns_boito = annual_eta(transversal_data=boito_fac_etas[0],
                          longitudinal_data=boito_fac_etas[1],
                          site=site, NS=True)
    ew_boito = annual_eta(transversal_data=boito_fac_etas[0],
                          longitudinal_data=boito_fac_etas[1],
                          site=site, NS=False)
    ##############################################################################################################

    # Specific reference #############################################################################################
    spec_ref_fac_etas = design_position_curvature(geom=geom,
                                                  curvature_design='SR').optical_analysis(cum_eff=cum_eff,
                                                                                          rec_aim=geom.rec_aim,
                                                                                          end_losses=False,
                                                                                          factorized=True)

    ns_spr = annual_eta(transversal_data=spec_ref_fac_etas[0],
                        longitudinal_data=spec_ref_fac_etas[1], site=site, NS=True)
    ew_spr = annual_eta(transversal_data=spec_ref_fac_etas[0],
                        longitudinal_data=spec_ref_fac_etas[1], site=site, NS=False)
    ##################################################################################################################

    # Non-uniform optimization routine ###############################################################################
    nun_or_file = Path(files_path, 'nun_or.json')

    if nun_or_file.is_file():
        with open(nun_or_file, 'r') as loaded_file:
            nun_or_dic = json.load(loaded_file)
            ns_nun_radii = array(nun_or_dic['ns'])
            ns_nun_lfr = lfr_design(geom=geom, radii=ns_nun_radii)
            ew_num_radii = array(nun_or_dic['ew'])
            ew_nun_lfr = lfr_design(geom=geom, radii=ew_num_radii)
    else:
        ns_nun_lfr, ew_nun_lfr = non_uniform_optimization_routine(geom=geom, site=site, cum_eff=cum_eff)
        ns_nun_radii = ns_nun_lfr.field.radius
        ew_nun_radii = ew_nun_lfr.field.radius
        nun_or_dic = {'ns': ns_nun_radii, 'ew': ew_nun_radii}
        dic2json(d=nun_or_dic, file_path=files_path, file_name='nun_or')

    ns_nun_fac_etas = ns_nun_lfr.optical_analysis(cum_eff=cum_eff)
    ns_nun = annual_eta(transversal_data=ns_nun_fac_etas[0],
                        longitudinal_data=ns_nun_fac_etas[1], site=site, NS=True)

    ew_nun_fac_etas = ew_nun_lfr.optical_analysis(cum_eff=cum_eff)
    ew_nun = annual_eta(transversal_data=ew_nun_fac_etas[0],
                        longitudinal_data=ew_nun_fac_etas[1], site=site, NS=False)
    ##################################################################################################################

    # output dictionaries ###############################
    ns_d = {'ns_dp': ns_dp,
            'ns_opt_dp': array([[ns_max_dp, ns_max_eta]]),
            'ns_boito': ns_boito,
            'ns_spr': ns_spr,
            'ns_sun': array([[ns_sun_dp, ns_sun_eta]]),
            'ns_nun': ns_nun}

    ew_d = {'ew_dp': ew_dp,
            'ew_opt_dp': array([[ew_max_dp, ew_max_eta]]),
            'ew_boito': ew_boito,
            'ew_spr': ew_spr,
            'ew_sun': array([[ew_sun_dp, ew_sun_eta]]),
            'ew_nun': ew_nun}
    ###################################################

    return ns_d, ew_d


def uniform_optimization_routine(geom: lfr_geometry, site: SiteData, source: RadialSource):
    files_path = Path(Path.cwd(), site.name, geom.name, source.name, 'efficiency')
    files_path.mkdir(parents=True, exist_ok=True)

    un_or_file = Path(files_path, 'un_or.json')

    if un_or_file.is_file():
        with open(un_or_file, 'r') as uploaded_file:
            un_dic = json.load(uploaded_file)

        ns_d, ew_d = un_dic['ns'], un_dic['ew']

        ns_d['ns_data'] = array(ns_d['ns_data'])
        ns_d['ns_opt'] = array(ns_d['ns_opt'])

        ew_d['ew_data'] = array(ew_d['ew_data'])
        ew_d['ew_opt'] = array(ew_d['ew_opt'])

    else:
        cum_eff = source.linear_cum_eff()

        fm_r = 2 * dst(geom.centers[0], geom.rec_aim)
        fm_l = 2 * dst(geom.centers[-1], geom.rec_aim)
        fm_max = max(fm_l, fm_r)

        radii = fm_max * linspace(start=0.5, stop=3., num=50).round(4)
        ns_etas = zeros(radii.shape[0])
        ew_etas = zeros(radii.shape[0])

        for i, r in enumerate(tqdm(radii)):
            fac_etas = lfr_design(geom=geom, radii=r).optical_analysis(cum_eff=cum_eff,
                                                                       factorized=True,
                                                                       end_losses=False,
                                                                       symmetric=False)

            ns_etas[i] = annual_eta(transversal_data=fac_etas[0], longitudinal_data=fac_etas[1], site=site, NS=True)
            ew_etas[i] = annual_eta(transversal_data=fac_etas[0], longitudinal_data=fac_etas[1], site=site, NS=False)

        ns = zeros(shape=(radii.shape[0], 2))
        ns.T[:] = radii, ns_etas

        ns_function = InterpolatedUnivariateSpline(*ns.T, k=4)
        ns_opt_radius = ns_function.derivative().roots()

        ew = zeros(shape=(radii.shape[0], 2))
        ew.T[:] = radii, ew_etas

        ew_function = InterpolatedUnivariateSpline(*ew.T, k=4)
        ew_opt_radius = ew_function.derivative().roots()

        ns_d = {'ns_data': ns,
                'ns_opt': array([ns_opt_radius[0], ns_function(ns_opt_radius)[0]])}

        ew_d = {'ew_data': ew,
                'ew_opt': array([ew_opt_radius[0], ew_function(ew_opt_radius)[0]])}

        out_dic = {'ns': ns_d, 'ew': ew_d}
        dic2json(d=out_dic, file_path=files_path, file_name='un_or')

    return ns_d, ew_d


def acceptance4designs(geom: lfr_geometry, site: SiteData, source: RadialSource):
    files_path = Path(Path.cwd(), site.name, geom.name, source.name, 'acceptance')
    files_path.mkdir(parents=True, exist_ok=True)

    acceptance_file = Path(files_path, 'acceptance_data.json')
    if acceptance_file.is_file():
        print(f'Acceptance data file already exist for geometry {geom.name} and source {source.name}.')
        print('Reading data file.')

        with open(acceptance_file, 'r') as loaded_file:
            acc_dic = json.load(loaded_file)

        ns_dic, ew_dic = acc_dic['ns'], acc_dic['ew']

        for k in ns_dic.keys():
            ns_dic[k] = array(ns_dic[k])

        for k in ew_dic.keys():
            ew_dic[k] = array(ew_dic[k])

    else:
        # Acceptance calculations data ####
        lv = 0.85
        dt = 0.05
        ref_value = 0.9
        ###################################

        # linear cumulative density function of the effective source ###
        cum_eff = source.linear_cum_eff()
        ################################################################

        # range of transversal_angles to be analyzed ################################
        theta_i, theta_f, dth = -85., 85., 2.5
        transversal_angles = arange(start=theta_i, stop=theta_f + dth, step=dth)
        #############################################################################

        print("Running Specific reference design acceptance calculations.")
        # Specific reference design ##################################################################################
        lfr_concentrator = design_position_curvature(geom=geom, curvature_design='SR')

        spr_acceptances = [lfr_concentrator.acceptance_angle(theta_t=theta, rec_aim=geom.rec_aim, cum_eff=cum_eff,
                                                             lv=lv, dt=dt, ref_value=ref_value)
                           for theta in tqdm(transversal_angles)]

        spr_data = zeros(shape=(transversal_angles.shape[0], 2))
        spr_data.T[:] = transversal_angles, spr_acceptances

        spr_eta = lfr_concentrator.optical_analysis(cum_eff=cum_eff, rec_aim=geom.rec_aim,
                                                    symmetric=False, factorized=True)[0].tolist()

        ###############################################################################################################

        # Sun reference design ########################################################################################
        ns, ew = non_uniform_design_analysis(geom=geom, site=site, source=source)
        ns_sun_dp = ns['ns_sun'].T[0][0]
        ew_sun_dp = ew['ew_sun'].T[0][0]

        print("Running NS Sun reference design acceptance calculations.")
        lfr_concentrator = design_position_curvature(geom=geom, curvature_design=ns_sun_dp)
        ns_sun_acceptances = [lfr_concentrator.acceptance_angle(theta_t=theta, rec_aim=geom.rec_aim, cum_eff=cum_eff,
                                                                lv=lv, dt=dt, ref_value=ref_value)
                              for theta in tqdm(transversal_angles)]

        ns_sun_data = zeros(shape=(transversal_angles.shape[0], 2))
        ns_sun_data.T[:] = transversal_angles, ns_sun_acceptances

        ns_sun_eta = lfr_concentrator.optical_analysis(cum_eff=cum_eff, rec_aim=geom.rec_aim,
                                                       symmetric=False, factorized=True)[0].tolist()

        print("Running EW Sun reference design acceptance calculations.")
        lfr_concentrator = design_position_curvature(geom=geom, curvature_design=ew_sun_dp)
        ew_sun_acceptances = [lfr_concentrator.acceptance_angle(theta_t=theta, rec_aim=geom.rec_aim, cum_eff=cum_eff,
                                                                lv=lv, dt=dt, ref_value=ref_value)
                              for theta in tqdm(transversal_angles)]

        ew_sun_data = zeros(shape=(transversal_angles.shape[0], 2))
        ew_sun_data.T[:] = transversal_angles, ew_sun_acceptances

        ew_sun_eta = lfr_concentrator.optical_analysis(cum_eff=cum_eff, rec_aim=geom.rec_aim,
                                                       symmetric=False, factorized=True)[0].tolist()

        print("Running Boito design acceptance calculations.")
        # Boito design ################################################################################################
        lfr_concentrator = boito_design(geom=geom, lat=deg2rad(site.lat))
        boito_acceptances = [lfr_concentrator.acceptance_angle(theta_t=theta, rec_aim=geom.rec_aim, cum_eff=cum_eff,
                                                               lv=lv, dt=dt, ref_value=ref_value)
                             for theta in tqdm(transversal_angles)]

        boito_data = zeros(shape=(transversal_angles.shape[0], 2))
        boito_data.T[:] = transversal_angles, boito_acceptances

        boito_eta = lfr_concentrator.optical_analysis(cum_eff=cum_eff, rec_aim=geom.rec_aim,
                                                      symmetric=False, factorized=True)[0].tolist()
        ###############################################################################################################

        # Non-uniform optimum #########################################################################################
        nun_or_file = Path(Path.cwd(), site.name, geom.name, source.name, 'efficiency', 'nun_or.json')
        with open(nun_or_file, 'r') as loaded_file:
            nun_or_dic = json.load(loaded_file)

        ns_nun_radii = array(nun_or_dic['ns'])
        ns_nun_lfr = lfr_design(geom=geom, radii=ns_nun_radii)

        print("Running NS NUN-OR design acceptance calculations.")
        ns_nun_acceptances = [ns_nun_lfr.acceptance_angle(theta_t=theta, rec_aim=geom.rec_aim, cum_eff=cum_eff,
                                                          lv=lv, dt=dt, ref_value=ref_value)
                              for theta in tqdm(transversal_angles)]

        ns_nun_data = zeros(shape=(transversal_angles.shape[0], 2))
        ns_nun_data.T[:] = transversal_angles, ns_nun_acceptances

        ns_nun_eta = ns_nun_lfr.optical_analysis(cum_eff=cum_eff, rec_aim=geom.rec_aim,
                                                 symmetric=False, factorized=True)[0].tolist()

        print("Running EW NUN-OR design acceptance calculations.")
        ew_nun_radii = array(nun_or_dic['ew'])
        ew_nun_lfr = lfr_design(geom=geom, radii=ew_nun_radii)

        ew_nun_acceptances = [ew_nun_lfr.acceptance_angle(theta_t=theta, rec_aim=geom.rec_aim, cum_eff=cum_eff,
                                                          lv=lv, dt=dt, ref_value=ref_value)
                              for theta in tqdm(transversal_angles)]

        ew_nun_data = zeros(shape=(transversal_angles.shape[0], 2))
        ew_nun_data.T[:] = transversal_angles, ew_nun_acceptances

        ew_nun_eta = ew_nun_lfr.optical_analysis(cum_eff=cum_eff, rec_aim=geom.rec_aim,
                                                 symmetric=False, factorized=True)[0].tolist()
        ###############################################################################################################

        # Uniform optimum #############################################################################################
        ns_un_dic, ew_un_dic = uniform_optimization_routine(geom=geom, site=site, source=source)

        print("Running NS UN-OR design acceptance calculations.")
        ns_un_lfr = lfr_design(geom=geom, radii=ns_un_dic['ns_opt'].T[0])
        ns_un_acceptances = [ns_un_lfr.acceptance_angle(theta_t=theta, rec_aim=geom.rec_aim, cum_eff=cum_eff,
                                                        lv=lv, dt=dt, ref_value=ref_value)
                             for theta in tqdm(transversal_angles)]

        ns_un_data = zeros(shape=(transversal_angles.shape[0], 2))
        ns_un_data.T[:] = transversal_angles, ns_un_acceptances

        ns_un_eta = ns_un_lfr.optical_analysis(cum_eff=cum_eff, rec_aim=geom.rec_aim,
                                               symmetric=False, factorized=True)[0].tolist()

        print("Running EW UN-OR design acceptance calculations.")
        ew_un_lfr = lfr_design(geom=geom, radii=ew_un_dic['ew_opt'].T[0])
        ew_un_acceptances = [ew_un_lfr.acceptance_angle(theta_t=theta, rec_aim=geom.rec_aim, cum_eff=cum_eff,
                                                        lv=lv, dt=dt, ref_value=ref_value)
                             for theta in tqdm(transversal_angles)]

        ew_un_data = zeros(shape=(transversal_angles.shape[0], 2))
        ew_un_data.T[:] = transversal_angles, ew_un_acceptances

        ew_un_eta = ew_un_lfr.optical_analysis(cum_eff=cum_eff, rec_aim=geom.rec_aim,
                                               symmetric=False, factorized=True)[0].tolist()

        ns_dic = {'boito': boito_data,
                  'boito_eta': boito_eta,

                  'spr': spr_data,
                  'spr_eta': spr_eta,

                  'sun': ns_sun_data,
                  'sun_eta': ns_sun_eta,

                  'nun': ns_nun_data,
                  'nun_eta': ns_nun_eta,

                  'un': ns_un_data,
                  'un_eta': ns_un_eta}

        ew_dic = {'boito': boito_data,
                  'boito_eta': boito_eta,

                  'spr': spr_data,
                  'spr_eta': spr_eta,

                  'sun': ew_sun_data,
                  'sun_eta': ew_sun_eta,

                  'nun': ew_nun_data,
                  'nun_eta': ew_nun_eta,

                  'un': ew_un_data,
                  'un_eta': ew_un_eta}

        acc_dic = {'ns': ns_dic, 'ew': ew_dic}
        dic2json(d=acc_dic, file_path=files_path, file_name='acceptance_data')

    return ns_dic, ew_dic


def averaged_acceptance(acceptance_data: array, eta_data: array, site: SiteData):

    # The acceptance angle as a function of the transversal incidence angle ###############
    acceptance_function = InterpolatedUnivariateSpline(*acceptance_data.T, k=3, ext=1)
    #######################################################################################

    # linear angles for particular location ####
    tt, tl = site.linear_angles()
    ############################################

    # Definition of average acceptance that does not account for the optical efficiency -- but only the DNI #########
    # ns_avg_acc = acceptance_function(tt).dot(site.tmy_data['dni'].values) / site.tmy_data['dni'].values.sum()
    # ew_avg_acc = acceptance_function(tl).dot(site.tmy_data['dni'].values) / site.tmy_data['dni'].values.sum()
    #################################################################################################################

    # Definition of average acceptance that does account for the optical efficiency #########
    eta_function = InterpolatedUnivariateSpline(*eta_data.T, k=1, ext=1)

    ns_flux = eta_function(tt) * site.tmy_data['dni'].values
    ns_avg_acc = acceptance_function(tt).dot(ns_flux) / ns_flux.sum()

    ew_flux = eta_function(tl) * site.tmy_data['dni'].values
    ew_avg_acc = acceptance_function(tl).dot(ew_flux) / ew_flux.sum()
    ########################################################################################

    return ns_avg_acc, ew_avg_acc


def avg_acc4designs(geom: lfr_geometry, site: SiteData, source: RadialSource):
    dic_keys = ['boito', 'spr', 'sun', 'nun', 'un']

    ns_acc, ew_acc = acceptance4designs(geom=geom, site=site, source=source)

    ns_dic = {}
    ew_dic = {}

    for k in dic_keys:
        ns_dic[k], ew_dic[k] = averaged_acceptance(acceptance_data=ns_acc[k], eta_data=ns_acc[f'{k}_eta'], site=site)

    return ns_dic, ew_dic


def check_geometry(geom: lfr_geometry, geom_file: Path):
    with open(geom_file, 'r') as loaded_file:
        loaded_geometry = json.load(loaded_file)

    signal = True if geom.export2dic() == loaded_geometry else False

    return signal


def check_source(source: RadialSource, source_file: Path):
    with open(source_file, 'r') as loaded_file:
        loaded_source = json.load(loaded_file)

    signal = True if source.export2dic() == loaded_source else False

    return signal


def curvature_designs_analysis(geom: lfr_geometry, site: SiteData, source: RadialSource):
    print(f'It is now {datetime.now().strftime("%H:%M:%S")}.')
    print(f'Running simulations for geometry {geom.name}, source {source.name}, in {site.name}.')

    geom_path = Path(Path.cwd(), f'{site.name}')
    geom_path.mkdir(parents=True, exist_ok=True)
    geom_file = Path(geom_path, f'{geom.name}_geometry.json')

    if geom_file.is_file():
        assert check_geometry(geom=geom, geom_file=geom_file), \
            'Geometry being analyzed does not match with the one already existing. Please, check it!'
    else:
        geom.export2json(file_path=geom_path)

    source_path = Path(geom_path, geom.name)
    source_path.mkdir(parents=True, exist_ok=True)
    source_file = Path(source_path, f'{source.name}.json')

    if source_file.is_file():
        assert check_source(source=source, source_file=source_file), \
            'Source being analyzed does not match with the one already existing. Please, check it!'
    else:
        source.export2json(file_path=source_path)

    non_uniform_data = non_uniform_design_analysis(geom=geom, site=site, source=source)
    uniform_data = uniform_optimization_routine(geom=geom, site=site, source=source)
    acceptance_data = acceptance4designs(geom=geom, site=site, source=source)
    avg_acceptance = avg_acc4designs(geom=geom, site=site, source=source)

    out_files_path = Path(source_path, source.name)
    ns, ew = non_uniform_data

    ns_etas = {'dp': ns['ns_dp'],
               'opt_dp': ns['ns_opt_dp'],
               'sun': ns['ns_sun'], 'spr': ns['ns_spr'],
               'boito': ns['ns_boito'],
               'nun': ns['ns_nun']}

    ew_etas = {'dp': ew['ew_dp'],
               'opt_dp': ew['ew_opt_dp'],
               'sun': ew['ew_sun'], 'spr': ew['ew_spr'],
               'boito': ew['ew_boito'],
               'nun': ew['ew_nun']}

    ns, ew = uniform_data
    ns_etas['un'] = ns['ns_opt']
    ns_etas['un_data'] = ns['ns_data']

    ew_etas['un'] = ew['ew_opt']
    ew_etas['un_data'] = ew['ew_data']

    out_dic = {'ns': ns_etas, 'ew': ew_etas}
    dic2json(d=out_dic, file_path=out_files_path, file_name='etas_data')

    return non_uniform_data, uniform_data, acceptance_data, avg_acceptance

########################################################################################################################
########################################################################################################################
