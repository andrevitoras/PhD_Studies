import random
from typing import Tuple

from numpy import array, absolute, pi, where, sqrt, ndarray, ones, zeros, log, power
from scipy.interpolate import interp1d

from solaren.niopy.geometric_transforms import dst, nrm
from solaren.scopy.linear_fresnel import Absorber, PrimaryField, PrimaryMirror, primaries_curvature_radius
from solaren.scopy.linear_fresnel import uniform_centers
from solaren.scopy.nio_concentrators import cpc_type, symmetric_cec2tube, symmetric_cec2evacuated_tube
from solaren.scopy.sunlight import SiteData, RadialSource

import seaborn
from matplotlib import pyplot as plt


########################################################################################################################
# Classes ##############################################################################################################


class VariableBounds:

    def __init__(self, height, width, gap, radius):
        self.height = height
        self.width = width
        self.gap = gap
        self.radius = radius


class EffectiveSource:

    def __init__(self,
                 name: str,
                 sunshape: str = None, size: float = None,
                 slope_error: float = 0., specular_error: float = 0.):
        # Auxiliary attributes
        self.name = name
        self.profile = sunshape
        self.sunshape_size = abs(size) if size is not None else 0.
        self.slope_error = abs(slope_error)
        self.specular_error = abs(specular_error)

        # overall surface error
        self.overall_error = sqrt(4 * self.slope_error ** 2 + self.specular_error ** 2)

        # RadialSource object to represent the sunshape profile
        self.source = RadialSource(profile=self.profile, size=self.sunshape_size)
        # cumulative effective source for the analytical computations
        self.cum_eff = self.source.linear_convolution(slope_error=self.slope_error,
                                                      specular_error=self.specular_error)


class BoundaryConditions:

    def __init__(self,
                 bounds: VariableBounds,
                 flux_threshold: float,
                 effective_source: EffectiveSource,
                 location: SiteData):
        self.bounds = bounds
        self.effective_source = effective_source
        self.location = location
        self.flux_threshold = abs(flux_threshold)


class lfc_geometry:

    def __init__(self,
                 receiver_height: float,
                 widths: ndarray,
                 centers: ndarray,
                 radii: ndarray,
                 force_sun_ref_radii: bool = False,
                 force_ipa_ref_radii: bool = False,
                 x: list = None, bounds: VariableBounds = None,
                 length=100000):

        assert widths.shape[0] == centers.shape[0], \
            'Number of mirrors does not stands for all variables! Please check it!'

        assert widths.shape[0] == radii.shape[0], \
            'Number of mirrors does not stands for all variables! Please check it!'

        self.decision_vector = array(x[1:]) if x is not None else None
        self.bounds = bounds

        # Evacuated absorber ##############################################################################
        # evacuated tube with a CEC secondary optic
        self.absorber = Absorber.evacuated_tube(center=array([0, receiver_height]),
                                                absorber_radius=35.0,
                                                outer_cover_radius=62.5,
                                                inner_cover_radius=62.5 - 3)

        # Mertins' (2009) receiver to test the cost function
        # mertins_tube_radius = 0.219 * 1e3 / 2
        # self.absorber = Absorber.evacuated_tube(center=array([0, receiver_height]),
        #                                         absorber_radius=mertins_tube_radius,
        #                                         outer_cover_radius=1.79*mertins_tube_radius,
        #                                         inner_cover_radius=1.79*mertins_tube_radius - 3)
        ##################################################################################################

        # CEC secondary optic design and flat absorber at CEC aperture ######################################
        self.primary_field_width = dst(centers[0], centers[-1]) + 0.5 * (widths[0] + widths[-1])

        self.secondary = self.add_cec_secondary()
        self.rec_aim = self.secondary.aperture_center
        self.s1 = self.secondary.contour[0]
        self.s2 = self.secondary.contour[-1]

        self.flat_absorber = Absorber.flat(width=dst(self.s1, self.s2),
                                           center=self.rec_aim,
                                           axis=nrm(self.s2 - self.s1))
        ######################################################################################################

        # Check if argument radii must be ignored to calculate curvature radius of each mirror
        # following sun reference design
        if force_sun_ref_radii:
            radiuses = primaries_curvature_radius(centers=centers,
                                                  rec_aim=self.rec_aim,
                                                  radius_design='zenithal')
        elif force_ipa_ref_radii:
            R_ipa = 2 * dst(centers[0], self.rec_aim)
            radiuses = ones(widths.shape[0]) * R_ipa
        else:
            radiuses = radii

        # The primary field
        self.nbr_pts = (widths * 300 / 750).max()
        heliostats = [PrimaryMirror(center=hc,
                                    width=w,
                                    radius=r,
                                    nbr_pts=int(self.nbr_pts))
                      for hc, w, r in zip(centers, widths, radiuses)]
        self.primary_field = PrimaryField(heliostats=heliostats)

        # Auxiliary attributes ######################################
        self.number_primaries = widths.shape[0]
        self.receiver_height = abs(receiver_height)
        self.length = abs(length)
        #############################################################

        # General attributes #######################################################################################
        self.net_area = self.primary_field.widths.sum() * self.length / 1.e6  # in m2
        self.gross_area = self.primary_field_width * self.length / 1.e6  # in m2
        self.absorber_area = 2 * pi * self.absorber.absorber_tube.radius * self.length / 1.e6  # in m2
        self.geometric_concentration = self.net_area / self.absorber_area
        self.dimensionless_parameters = self.dimensionless_parameters()
        ############################################################################################################

    def is_uniform_width(self):
        return True if self.primary_field.widths.std().round(5) == 0. else False

    def is_uniform_radius(self):
        return True if self.primary_field.radius.std().round(5) == 0. else False

    def is_uniform_shift(self):
        return True if self.primary_field.shifts.std().round(5) == 0. else False

    def dimensionless_parameters(self):

        filling_factor = self.net_area / self.gross_area
        aspect_factor = self.primary_field_width / (2 * self.receiver_height)

        if self.is_uniform_width() and self.is_uniform_shift():
            params = filling_factor, aspect_factor, self.primary_field.shifts.mean() / self.primary_field.widths.mean()
        else:
            params = filling_factor, aspect_factor, 0.

        return params

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

    def optical_analysis(self,
                         cum_eff: array,
                         end_losses: bool = True,
                         symmetric: bool = True):

        transversal_data, longitudinal_data = self.primary_field.optical_analysis(flat_absorber=self.flat_absorber,
                                                                                  length=self.length,
                                                                                  rec_aim=self.rec_aim,
                                                                                  cum_eff=cum_eff,
                                                                                  end_losses=end_losses,
                                                                                  symmetric=symmetric,
                                                                                  factorized=True)

        return transversal_data, longitudinal_data

    def energy_collection_factor(self,
                                 cum_eff: array,
                                 location: SiteData,
                                 orientation='NS',
                                 end_losses: bool = True,
                                 symmetric: bool = True,
                                 flux_threshold: float = 0.):

        transversal_data, longitudinal_data = self.optical_analysis(cum_eff=cum_eff,
                                                                    end_losses=end_losses,
                                                                    symmetric=symmetric)

        # factorized efficiency functions -- functions here are vectorized.
        eta_0 = longitudinal_data[0][-1]
        transversal_function = interp1d(*transversal_data.T, kind='linear')
        longitudinal_function = interp1d(*longitudinal_data.T, kind='linear')

        # transversal and longitudinal-solar incidence angles
        theta_t, _, theta_ls = location.linear_angles(NS=True if orientation == 'NS' else False,
                                                      solar_longitudinal=True)
        if symmetric:
            location.tmy_data['theta_t'] = absolute(theta_t)
        else:
            location.tmy_data['theta_t'] = theta_t
        location.tmy_data['theta_ls'] = absolute(theta_ls)

        # selecting only the hourly data with dni > 0
        df = location.tmy_data[location.tmy_data['dni'] > 0].copy(deep=True)

        # optical efficiency for each hour
        df['efficiency'] = transversal_function(df['theta_t']) * longitudinal_function(df['theta_ls']) / eta_0
        # the flux at the absorber for each hour
        df['absorber_flux'] = df['efficiency'] * df['dni'] * self.net_area

        # the useful flux is only for values above a particular threshold
        df['useful_flux'] = df['absorber_flux'] - 1000 * flux_threshold * self.absorber_area
        df['useful_flux'] = where(df['useful_flux'] > 0, df['useful_flux'], 0)

        # the energy collection factor -- a useful annual efficiency
        # sum of useful flux to the sum of available flux (dni * aperture area) ratio
        ecf = df['useful_flux'].sum() / (df['dni'].sum() * self.net_area)

        return ecf

    def specific_cost(self):

        cost = mertins_specific_direct_cost(
            widths=self.primary_field.widths / 1000.,
            gaps=self.primary_field.gaps / 1000.,
            tube_diameter=(2 * self.absorber.absorber_tube.radius) / 1000.,
            receiver_height=self.receiver_height / 1000.,
            dh=4.0)

        return cost


########################################################################################################################
########################################################################################################################


########################################################################################################################
# Auxiliary functions ##################################################################################################

"""
This first set of functions implement Mertins' [1] model to estimate the direct specific cost of a linear Fresnel solar
field.

[1] Mertins, M., 2009. Technische und wirtschaftliche Analyse von horizontalen Fresnel-Kollektoren. 
    University of Karlsruhe, PhD Thesis.
"""


def mertins_specific_direct_cost(
        widths: array,
        gaps: array,
        tube_diameter: float,
        receiver_height: float,
        dh=4.0):

    # number of primary mirrors in the field
    n = widths.shape[0]

    assert n == gaps.shape[0] + 1, "The number of mirrors and gaps does not fit. Please, verify!"

    Hr = abs(receiver_height)  # receiver height
    da = abs(tube_diameter)

    Cm = array([mirror_cost_factor(mirror_width=w) for w in widths])
    Cg = array([gap_cost_factor(mirror_gap=g) for g in gaps])

    Ce = elevation_cost_factor(tube_diameter=da)
    Cr = receiver_cost_factor(tube_diameter=da)

    c = (Cm.sum() + Cg.dot(gaps) + Ce * (Hr + dh) + Cr) / widths.sum()

    return c


def mirror_cost_factor(mirror_width: float,
                       cm0=30.5,
                       w0=0.5,
                       nm=1.) -> float:

    w = abs(mirror_width)
    cm = cm0 * (w/w0)**nm

    return cm


def gap_cost_factor(mirror_gap: float,
                    cg0=11.5,
                    g0=0.01,
                    ng=1.) -> float:

    g = abs(mirror_gap)
    cg = cg0 * (g / g0) ** ng

    return cg


def mean_log_cost_factor(da: float,
                         d0: float,
                         c0: float,
                         components_costs: array,
                         components_scale_factors: array) -> float:
    """

    :param da:
    :param d0:
    :param c0:
    :param components_scale_factors:
    :param components_costs:
    :return:
    """

    # num = [(ci / c0) * (d / d0)**ni
    #        for ci, ni in zip(components_costs, components_scale_factors)]

    num = log((components_costs / c0).dot(power(da/d0, components_scale_factors)))
    den = 1 if da == d0 else log(da / d0)

    sf = num / den
    c = c0 * (da / d0) ** sf

    return c


def elevation_cost_factor(tube_diameter: float,
                          d0=0.219,
                          ce0=19.8,
                          n_ci_values=(1.4, 1., 1.),
                          ci_values=(14.2, 0.9, 4.6)) -> float:

    ce = mean_log_cost_factor(da=tube_diameter,
                              d0=d0,
                              c0=ce0,
                              components_costs=array(ci_values),
                              components_scale_factors=array(n_ci_values))

    return ce


def receiver_cost_factor(tube_diameter: float,
                         d0=0.219,
                         cr0=654.0,
                         n_ci_values=(2, 0.9, 0.7, 1.4, 0.6, 0.6),
                         ci_values=(161.2, 56.6, 116.4, 136.5, 26.4, 112.6)) -> float:

    cr = mean_log_cost_factor(da=tube_diameter,
                              d0=d0,
                              c0=cr0,
                              components_costs=array(ci_values),
                              components_scale_factors=array(n_ci_values))

    return cr

########################################################################################################################


def normalize(a: float, a0: float, a1: float):
    """Normalizes the value 'a' within the range [a0, a1] to [0, 1]."""
    return (a - a0) / (a1 - a0)


def denormalize(x: float, a0: float, a1: float) -> float:
    """Converts the normalized value 'x' from [0, 1] back to the range [a0, a1]."""
    return x * (a1 - a0) + a0


def number_of_decision_variables(number_of_mirrors: int, configuration='non-uniform'):

    if configuration == 'non-uniform':
        nbr_mirrors_variables = number_of_mirrors // 2 if number_of_mirrors % 2 == 0 else 1 + (number_of_mirrors // 2)
        nbr_dec_variables = 3 * nbr_mirrors_variables + 1
    elif configuration == 'uniform':
        nbr_dec_variables = 4
    elif configuration == 'variable-radius':
        nbr_mirrors_variables = number_of_mirrors // 2 if number_of_mirrors % 2 == 0 else 1 + (number_of_mirrors // 2)
        nbr_dec_variables = 3 + nbr_mirrors_variables
    elif configuration == 'nun-sun-ref':
        nbr_mirrors_variables = number_of_mirrors // 2 if number_of_mirrors % 2 == 0 else 1 + (number_of_mirrors // 2)
        nbr_dec_variables = 2 * nbr_mirrors_variables + 1
    elif configuration == 'un-ipa-ref':
        nbr_dec_variables = 3
    else:
        raise ValueError('Please, input a value primary field configuration!')

    return nbr_dec_variables


def normalized_non_uniform_decision_vector(number_of_mirrors: int) -> list:
    # Calculating the number of decision variables as function of the number of primary mirrors
    n_decision_variables = number_of_decision_variables(number_of_mirrors=number_of_mirrors)

    # generating the decision vector with the normalized variables
    x = [number_of_mirrors] + [random.uniform(a=0, b=1) for i in range(n_decision_variables)]

    return x


def normalized_uniform_decision_vector(number_of_mirrors: int) -> list:
    # Calculating the number of decision variables as function of the number of primary mirrors
    n_decision_variables = number_of_decision_variables(number_of_mirrors=number_of_mirrors,
                                                        configuration='uniform')

    # generating the decision vector with the normalized variables
    x = [number_of_mirrors] + [random.uniform(a=0, b=1) for i in range(n_decision_variables)]

    return x


def normalized_variable_radius_decision_vector(number_of_mirrors: int) -> list:
    # Calculating the number of decision variables as function of the number of primary mirrors
    n_decision_variables = number_of_decision_variables(number_of_mirrors=number_of_mirrors,
                                                        configuration='variable-radius')

    # generating the decision vector with the normalized variables
    x = [number_of_mirrors] + [random.uniform(a=0, b=1) for i in range(n_decision_variables)]

    return x


def normalized_nun_sun_ref_decision_vector(number_of_mirrors: int) -> list:
    # Calculating the number of decision variables as function of the number of primary mirrors
    n_decision_variables = number_of_decision_variables(number_of_mirrors=number_of_mirrors,
                                                        configuration='nun-sun-ref')

    # generating the decision vector with the normalized variables
    x = [number_of_mirrors] + [random.uniform(a=0, b=1) for i in range(n_decision_variables)]

    return x


def normalized_un_ipa_ref_decision_vector(number_of_mirrors: int) -> list:
    # Calculating the number of decision variables as function of the number of primary mirrors
    n_decision_variables = number_of_decision_variables(number_of_mirrors=number_of_mirrors,
                                                        configuration='un-ipa-ref')

    # generating the decision vector with the normalized variables
    x = [number_of_mirrors] + [random.uniform(a=0, b=1) for i in range(n_decision_variables)]

    return x


def denormalized_non_uniform_variables(x: list,
                                       bounds: VariableBounds):
    number_of_mirrors = x[0]

    n = (len(x) - 2) // 3  # number of decision variables

    x1 = x[1]  # the normalized values of the receiver height
    x2 = x[2: n + 2]  # normalized vector of the mirrors widths
    x3 = x[n + 2: 2 * n + 2]  # normalized vector mirrors radii
    x4 = x[2 * n + 2:]  # normalized vector mirrors gaps

    receiver_height = denormalize(x=x1, a0=bounds.height[0], a1=bounds.height[1])
    widths = [denormalize(x=v, a0=bounds.width[0], a1=bounds.width[1]) for v in x2]
    radii = [denormalize(x=v, a0=bounds.radius[0], a1=bounds.radius[1]) for v in x3]
    gaps = [denormalize(x=v, a0=bounds.gap[0], a1=bounds.gap[1]) for v in x4]

    return number_of_mirrors, receiver_height, widths, radii, gaps


def non_uniform_geometry(x: list,
                         bounds: VariableBounds) -> lfc_geometry:

    number_of_mirrors, receiver_height, widths, radii, gaps = denormalized_non_uniform_variables(x=x, bounds=bounds)

    if number_of_mirrors % 2 != 0:
        centers = [[0, 0]]
    else:
        centers = [[gaps[0] + 0.5 * widths[0], 0]]

    for i in range(len(widths) - 1):
        centers.append([centers[i][0] + gaps[i + 1] + 0.5 * (widths[i] + widths[i + 1]), 0])

    centers = array(centers)

    if number_of_mirrors % 2 != 0:
        widths = widths[::-1] + widths[1:]
        radii = radii[::-1] + radii[1:]
        centers = centers[::-1].tolist() + (-centers[1:]).tolist()
    else:
        widths = widths[::-1] + widths
        radii = radii[::-1] + radii
        centers = centers[::-1].tolist() + (-centers).tolist()

    widths = array(widths)
    radii = array(radii)
    centers = array(centers)

    # the geometry object
    solution = lfc_geometry(x=x,
                            bounds=bounds,
                            widths=widths,
                            centers=centers,
                            radii=radii,
                            receiver_height=receiver_height)

    return solution


def denormalized_uniform_variables(x: list,
                                   bounds: VariableBounds):
    number_of_mirrors = x[0]

    x1 = x[1]  # the normalized values of the receiver height
    x2 = x[2]  # normalized vector of the mirrors width
    x3 = x[3]  # normalized vector mirrors radius
    x4 = x[4]  # normalized vector mirrors gap

    receiver_height = denormalize(x=x1, a0=bounds.height[0], a1=bounds.height[1])
    width = denormalize(x=x2, a0=bounds.width[0], a1=bounds.width[1])
    radius = denormalize(x=x3, a0=bounds.radius[0], a1=bounds.radius[1])
    gap = denormalize(x=x4, a0=bounds.gap[0], a1=bounds.gap[1])

    return number_of_mirrors, receiver_height, width, radius, gap


def uniform_geometry(x: list, bounds: VariableBounds):
    number_of_mirrors, receiver_height, width, radius, gap = denormalized_uniform_variables(x=x, bounds=bounds)

    # the distance between two centers
    s = gap + width

    # the center points of the uniform configuration
    centers = uniform_centers(mirror_width=width, center_distance=s, nbr_mirrors=number_of_mirrors)

    widths = ones(number_of_mirrors) * width
    radii = ones(number_of_mirrors) * radius

    solution = lfc_geometry(x=x,
                            bounds=bounds,
                            widths=widths,
                            centers=centers,
                            radii=radii,
                            receiver_height=receiver_height)

    return solution


def denormalized_var_radius_variables(x: list,
                                      bounds: VariableBounds) -> Tuple[int, float, float, float, list]:
    number_of_mirrors = x[0]
    x1 = x[1]  # the normalized values of the receiver height
    x2 = x[2]  # normalized vector of the mirrors width
    x3 = x[3]  # normalized vector mirrors gap
    x4 = x[4:]  # normalized vector mirrors radii

    receiver_height = denormalize(x=x1, a0=bounds.height[0], a1=bounds.height[1])
    width = denormalize(x=x2, a0=bounds.width[0], a1=bounds.width[1])
    gap = denormalize(x=x3, a0=bounds.gap[0], a1=bounds.gap[1])
    radii = [denormalize(x=x, a0=bounds.radius[0], a1=bounds.radius[1])
             for x in x4]

    return number_of_mirrors, receiver_height, width, gap, radii


def variable_radius_geometry(x: list, bounds: VariableBounds):
    number_of_mirrors, receiver_height, width, gap, radii = denormalized_var_radius_variables(x=x, bounds=bounds)

    # the distance between two centers
    s = gap + width

    # the center points of the uniform configuration
    centers = uniform_centers(mirror_width=width, center_distance=s, nbr_mirrors=number_of_mirrors)
    widths = ones(number_of_mirrors) * width

    if number_of_mirrors % 2 != 0:
        radiuses = radii[::-1] + radii[1:]
    else:
        radiuses = radii[::-1] + radii

    radiuses = array(radiuses)

    solution = lfc_geometry(x=x,
                            bounds=bounds,
                            widths=widths,
                            centers=centers,
                            radii=radiuses,
                            receiver_height=receiver_height)

    return solution


def denormalized_nun_sun_ref_variables(x: list,
                                       bounds: VariableBounds):
    number_of_mirrors = x[0]
    n = (len(x) - 2) // 2  # number of decision variables

    x1 = x[1]  # normalized value of the receiver height
    x2 = x[2: n + 2]  # normalized vector of the mirrors widths
    x3 = x[n + 2:]  # normalized vector mirrors gaps

    receiver_height = denormalize(x=x1, a0=bounds.height[0], a1=bounds.height[1])
    widths = [denormalize(x=v, a0=bounds.width[0], a1=bounds.width[1]) for v in x2]
    gaps = [denormalize(x=v, a0=bounds.gap[0], a1=bounds.gap[1]) for v in x3]

    return number_of_mirrors, receiver_height, widths, gaps


def nun_sun_ref_geometry(x: list,
                         bounds: VariableBounds) -> lfc_geometry:

    number_of_mirrors, receiver_height, widths, gaps = denormalized_nun_sun_ref_variables(x=x, bounds=bounds)

    if number_of_mirrors % 2 != 0:
        centers = [[0, 0]]
    else:
        centers = [[gaps[0] + 0.5 * widths[0], 0]]

    for i in range(len(widths) - 1):
        centers.append([centers[i][0] + gaps[i + 1] + 0.5 * (widths[i] + widths[i + 1]), 0])

    centers = array(centers)

    if number_of_mirrors % 2 != 0:
        widths = widths[::-1] + widths[1:]
        centers = centers[::-1].tolist() + (-centers[1:]).tolist()
    else:
        widths = widths[::-1] + widths
        centers = centers[::-1].tolist() + (-centers).tolist()

    widths = array(widths)
    centers = array(centers)
    radii = zeros(widths.shape[0])

    # the geometry object
    solution = lfc_geometry(x=x,
                            bounds=bounds,
                            receiver_height=receiver_height,
                            widths=widths,
                            centers=centers,
                            radii=radii,
                            force_sun_ref_radii=True)

    return solution


def denormalized_un_ipa_ref_variables(x: list,
                                      bounds: VariableBounds):
    number_of_mirrors = x[0]

    x1 = x[1]  # normalized value of the receiver height
    x2 = x[2]  # normalized mirrors width
    x3 = x[3]  # normalized mirrors gap

    receiver_height = denormalize(x=x1, a0=bounds.height[0], a1=bounds.height[1])
    width = denormalize(x=x2, a0=bounds.width[0], a1=bounds.width[1])
    gap = denormalize(x=x3, a0=bounds.gap[0], a1=bounds.gap[1])

    return number_of_mirrors, receiver_height, width, gap


def un_ipa_ref_geometry(x: list,
                        bounds: VariableBounds) -> lfc_geometry:

    number_of_mirrors, receiver_height, width, gap = denormalized_un_ipa_ref_variables(x=x, bounds=bounds)

    # the distance between two centers
    s = gap + width

    # the center points of the uniform configuration
    centers = uniform_centers(mirror_width=width, center_distance=s, nbr_mirrors=number_of_mirrors)

    widths = ones(number_of_mirrors) * width
    radii = zeros(widths.shape[0])

    solution = lfc_geometry(x=x,
                            bounds=bounds,
                            widths=widths,
                            centers=centers,
                            radii=radii,
                            receiver_height=receiver_height,
                            force_ipa_ref_radii=True)

    return solution


def cec4tube(primary_field_width: float,
             tube_center: array,
             tube_radius: float,
             gap_radius: float = None,
             points_per_section: int = 50):
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


def SchottTubeHeatLosses(Tabs: float):
    """
    This function compute the heat losses per unit of linear meter (W/m) of the 70 mm diameter SCHOTT evacuated tube.
    It uses the correlation function from Burkholder and Kutscher [1] that relates absorber temperature and heat losses.

    [1] Burkholder, F., Kutscher, C., 2009. Heat Loss Testing of Schott’s 2008 PTR70 Parabolic Trough Receiver.
    Golden, CO (United States). https://doi.org/10.2172/1369635.

    :param Tabs: absorber tube temperature, in ºC

    :return: the heat losses in W/m of linear tube.
    """

    return 0.141*Tabs + 6.48e-9 * (Tabs**4)

########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################


if __name__ == '__main__':
    seaborn.set(style='whitegrid')
    seaborn.set_context('notebook')

    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'NewComputerModern10'

    ticks_font_size = 14
    plt.rcParams['xtick.labelsize'] = ticks_font_size
    plt.rcParams['ytick.labelsize'] = ticks_font_size

    axes_label_size = 14
    plt.rcParams['axes.labelsize'] = axes_label_size
    plt.rcParams['axes.titlesize'] = axes_label_size

    #######################################################################

    var_bounds = VariableBounds(height=(4 * 1e3, 20 * 1e3),
                                width=(0.2 * 1e3, 2 * 1e3),
                                gap=(0, 2 * 1e3),
                                radius=(0, 100 * 1e3))

    dec_vec = normalized_un_ipa_ref_decision_vector(number_of_mirrors=25)
    lfc = un_ipa_ref_geometry(x=dec_vec, bounds=var_bounds)

    # d = normalized_uniform_decision_vector(number_of_mirrors=12)
    # lfc = uniform_geometry(x=d, bounds=var_bounds)

    # Mertins (2009) geometry ###########################################################################
    # width = 0.5 * 1000
    # gap = 0.01 * 1000
    # centers = uniform_centers(mirror_width=width, center_distance=gap + width, nbr_mirrors=48)
    # widths = ones(48) * width
    # lfc = lfc_geometry(widths=widths, centers=centers, receiver_height=9.0 * 1000,
    #                    radii=zeros(48))
    #
    ########################################################################################################

    fig = plt.figure(dpi=300)
    ax = fig.add_subplot()

    [ax.plot(*hel.contour.T / 1000,
             color='black', lw=1.5, label='Primary mirrors' if i == 0 else None)
     for i, hel in enumerate(lfc.primary_field.heliostats)]

    [ax.plot(*hel.center.T / 1000,
             color='orange', lw=0, marker='.', ms=4, label='Centers' if i == 0 else None)
     for i, hel in enumerate(lfc.primary_field.heliostats)]

    ax.plot(*lfc.absorber.absorber_tube.contour.T / 1000, color='red', lw=0.5, label='Absorber tube')
    ax.plot(*lfc.secondary.contour.T / 1000, color='blue', label='CEC secondary', lw=0.5)

    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$z$ [m]')
    ax.set_title('Linear Fresnel Collector: general view')

    ax.legend()
    plt.axis('equal')
    plt.tight_layout()
    plt.show()

    # fig = plt.figure(dpi=300)
    # ax = fig.add_subplot()
    #
    # ax.plot(*lfc.absorber.absorber_tube.contour.T / 1000, color='red', lw=0.75, label='Absorber tube')
    # ax.plot(*lfc.absorber.outer_tube.contour.T / 1000, color='green', lw=0.75, label='Glass cover')
    # ax.plot(*lfc.absorber.inner_tube.contour.T / 1000, color='green', lw=0.75)
    #
    # ax.plot(*lfc.secondary.contour.T / 1000, color='blue', label='CEC secondary')
    #
    # ax.set_xlabel('$x$ [m]')
    # ax.set_ylabel('$z$ [m]')
    # ax.set_title('Linear Fresnel Collector: receiver view')
    #
    # ax.legend()
    # plt.axis('equal')
    # plt.tight_layout()
    # plt.show()
