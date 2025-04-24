"""
created by André Santos (andrevitoras@gmail.com / avas@uevora.pt)

The set of functions presented here considers a transversal plane analysis of linear Fresnel collectors.
It is based on a XY plane, so that point and vectors lie on such a plane.

Thus, a positive transversal incidence angle refers to positive rotation regarding the Z-axis.
In the same sense, a positive angular position also refers to a positive rotation regarding the Z-axis.

Here, the term heliostat refers to a primary mirror of a linear Fresnel collector.

"""
import pandas as pd
from numpy import linspace, array, tan, arange, ones, zeros, deg2rad, power, pi, arctan, cos, sign, sqrt
from scipy.optimize import fsolve

from solaren.niopy.geometric_transforms import V, nrm, dst, R, ang
from solaren.niopy.plane_curves import PlaneCurve, par
from utils import rmse


########################################################################################################################
# Classes and functions for this program ###############################################################################

def cylindrical_curvature(center: array, aim: array, theta_d: float = 0.0) -> float:
    """
    A function to calculate the ideal curvature radius of a cylindrical heliostat as defined by Rabl [1, p.179].

    :param center: heliostat's center point.
    :param aim: aim point at the receiver.
    :param theta_d: design position, a transversal incidence angle (in degrees).
    :return: This function returns the ideal cylindrical curvature.

    References:
    [1] Rabl A. Active Solar Collectors and Their Applications. New York: Oxford University Press; 1985.

    It is important to highlight that calculations are for a xy-plane, where transversal incidence angles are positive
    on the left side of the y-axis direction (a positive rotation about the z-axis). The same definition is used
    for the heliostat angular position computed by the relation between the center point and the aim-point at the
    receiver.

    """

    # Angle from the horizontal which defines the direction of the incoming sunlight at the transversal plane.
    alpha = 0.5 * pi + deg2rad(theta_d)
    # the incidence vector in the transversal plane
    vi = V(alpha)

    # Forcing the center and aim as 2D array points: [x, y]
    hc = array([center[0], center[-1]])
    f = array([aim[0], aim[-1]])

    # The focusing vector
    vf = f - hc

    # It is better to check for a specific reference design and calculating it directly.
    # Check if the direction of the incoming sunlight is aligned with the mirror focusing vector since
    # the function 'ang(u, v)' used here sometimes calculates a wrong value when u || v.
    # Then, calculate the curvature absorber_radius.
    if abs(vf.dot(vi).round(5)) == 1:
        r = 2 * dst(hc, f)
    else:
        mi = 0.5 * ang(f - hc, vi)
        r = 2. * dst(hc, f) / cos(mi)

    return r


def angular_position(center: array, aim: array) -> float:
    """
    This function calculate the heliostat angular position

    :param center:
    :param aim:
    :return:
    """

    Iy = array([0, 1])
    vf = aim - center

    lamb = sign(center[0]) * ang(vf, Iy)

    return lamb


def tracking_angle(center: array, aim: array, theta_t_rad: float) -> float:

    lamb = angular_position(center=center, aim=aim)
    tau = 0.5 * (theta_t_rad + lamb)

    return tau


def design_flat_mirror(hc: array, w: float, nbr_pts: int):
    """
    This function returns the surface points of a flat shape heliostat.
    This set of points define the heliostat contour.

    Units should be in millimeters to be coherent with the classes, methods, and functions presented in this module.

    :param hc: heliostat center point
    :param w: heliostat width
    :param nbr_pts: number of point to parametrize.

    :return: This function returns a list of points from the function of the heliostat.
    """

    n_pts = nbr_pts + 1 if nbr_pts % 2 == 0 else nbr_pts
    center = array([hc[0], hc[-1]])

    hel_pts = zeros(shape=(n_pts, 2))
    hel_pts[:, 0] = linspace(start=-0.5 * w, stop=0.5 * w, num=n_pts) + center[0]

    return hel_pts


def design_cylindrical_mirror(center: array, width: float, radius: float, nbr_pts: int):

    # The array of x-values which the heliostat ranges.
    x_range = linspace(start=-0.5 * width, stop=+0.5 * width, num=nbr_pts)

    # Ensure that the center point is a XY array point.
    hc = array([center[0], center[-1]])

    # the function which analytically describes the cylindrical surface which comprises the heliostat
    def y(x): return -sqrt(radius ** 2 - x ** 2) + radius

    # the computation of the points which discretize the heliostat surface
    hel_pts = array([[x, y(x)] for x in x_range]) + hc

    return hel_pts


def design_nested_parabolic_mirror(center: array, focal_length: float, width: float, n_pts: int):

    nbr_pts = n_pts if n_pts % 2 != 0 else n_pts + 1

    hc = array([center[0], center[-1]])
    w = abs(width)
    fm = abs(focal_length)

    hel_pts = zeros(shape=(nbr_pts, 2))
    x_range = linspace(start=-0.5 * w, stop=0.5 * w, num=nbr_pts)

    def par_f(x): return power(x, 2) / (4 * fm)

    hel_pts[:] = [[v, par_f(v)] for v in x_range]

    return hel_pts + hc


def design_parabolic_mirror(center: array, aim: array, width: float, nbr_pts: int, theta_d=0.):
    """
    This function (...)

    :param center:
    :param aim:
    :param width:
    :param nbr_pts:
    :param theta_d: The design position, in degrees

    :return: A set of [x,y] points that defined the contour of the mirror
    """

    # design position in radians
    theta_d_rad = deg2rad(theta_d)

    # The number of points to discretize the parabolic heliostat contour.
    # It must be an odd number for the heliostat center be a point in the array
    n_pts = nbr_pts + 1 if nbr_pts % 2 == 0 else nbr_pts

    # Picking values and changing vectors to a 2D dimension
    w = abs(width)
    hc = array([center[0], center[-1]])
    f = array([aim[0], aim[-1]])

    # calculates the tracking angle
    tau = tracking_angle(center=center, aim=aim, theta_t_rad=theta_d_rad)

    # the angle which the parabola optical axis makes with the horizontal axis
    alpha = 0.5 * pi + theta_d_rad
    optical_axis = V(alpha)

    # Check for the specific reference design for a nested parabola design or not
    if nrm(f - hc).dot(optical_axis).round(3) == 1.0:
        hel_pts = design_nested_parabolic_mirror(center=hc, focal_length=dst(f, hc), width=w, n_pts=n_pts)
    else:
        # parabola function to compute the vertex point
        par_f = par(alpha=alpha,
                    f=f,
                    p=hc)

        vertex = par_f(pi)

        # nested parabola focal distance. For the calculation in the rotated reference frame.
        fm = dst(vertex, f)

        # nested parabola equation.
        def centered_par(x):
            return (x ** 2) / (4 * fm)

        # Map the heliostat center in the rotated reference frame
        # In this sense, a tilted parabola becomes a nested one
        pn = R(alpha=-theta_d_rad).dot(hc - vertex)

        # Derivative of the nested parabola for the mapped point in the rotated frame.
        d_pn = pn[0] / (2 * fm)

        # points distant half-width of the mapped point in the rotated frame.
        loc_slope = arctan(d_pn)
        p1 = pn + 0.5 * w * V(loc_slope)
        p2 = pn - 0.5 * w * V(loc_slope)

        # The edge points which define the parabolic heliostat are calculated by the interception between
        # the nested parabola equation and straight lines which are normal to the tangent at the
        # heliostat center mapped at the rotated reference frame, i.e., point 'pn'.

        def sl1(z):
            return p1[1] + (z - p1[0]) * tan(0.5 * pi + loc_slope)

        def sl2(z):
            return p2[1] + (z - p2[0]) * tan(0.5 * pi + loc_slope)

        e1 = fsolve(lambda z: centered_par(z) - sl1(z), pn[0])[0]
        e2 = fsolve(lambda z: centered_par(z) - sl2(z), pn[0])[0]

        # creates the range of x values comprised between the x components of the edges, i.e., 'e1' and 'e2'
        x_values = linspace(start=min(e1, e2), stop=max(e1, e2), num=n_pts)

        # Calculate the array of points which comprises the heliostat in the rotated frame
        rot_pts = zeros(shape=(n_pts, 2))
        rot_pts[:] = [[x, centered_par(x)] for x in x_values]

        # Calculating the points in the fixed reference frame
        rm = R(alpha=theta_d_rad)
        rot_pts = rm.dot(rot_pts.T).T + vertex

        # Rotating the heliostat to the horizontal position
        rm = R(alpha=-tau)
        hel_pts = rm.dot((rot_pts - hc).T).T + hc

    return hel_pts


class ParabolicHeliostat:
    """
    This class aim to represent primary mirrors with a parabolic shape. There are two design possibilities:
    (1) A design was proposed by Häberle [1], which considers a design position and an aim point;
    (2) the vertical design, which considers just the focal length of the parabola, here define by the boolean argument
    'forced_design'.


    [1] Häberle A. Linear Fresnel Collectors. Solar Energy, New York, NY: Springer New York; 2013, p. 72–8.
    https://doi.org/10.1007/978-1-4614-5806-7_679.

    """

    def __init__(self, center: array, width: float, theta_d: float = None, aim: array = None,
                 forced_design=False, focal_length: float = None,
                 nbr_pts=121):
        """

        :param center:
        :param width:
        :param theta_d: The design position, in degrees.
        :param aim:
        :param forced_design:
        :param focal_length:
        :param nbr_pts:
        """

        self.shape = 'parabolic'
        self.width = abs(width)

        self.center = array([center[0], center[-1]])

        # Ensure an odd number of points
        self.n_pts = nbr_pts + 1 if nbr_pts % 2 == 0 else nbr_pts

        if not forced_design and theta_d is not None and aim is not None:
            aim_pt = array([aim[0], aim[-1]])
            self.contour = design_parabolic_mirror(center=self.center, width=self.width, nbr_pts=self.n_pts,
                                                   aim=aim_pt, theta_d=theta_d)

        elif forced_design and focal_length is not None:
            self.contour = design_nested_parabolic_mirror(center=self.center, focal_length=focal_length,
                                                          width=self.width, n_pts=self.n_pts)

        else:
            raise ValueError('Invalid arguments. Please, see ParabolicHeliostat class documentation.')

        # Attribute that holds the Heliostat object as a PlaneCurve object.
        self.curve = self.as_plane_curve()
        self.normals = self.curve.normals2surface()

        self.seg_pts, self.seg_normals = self.segments_of_equal_projected_aperture()

    def segments_of_equal_projected_aperture(self):
        """
        This method calculates [x,0, z] point and normal vectors in the heliostat contour that are the central points
        of segments that has the same projected width in the aperture of the heliostat.
        These points ensure a uniform discretization of the mirror aperture and are essential to compute efficiency
        calculations.

        :return: A tuple of [x,0,z] point-arrays and [x,0,z] vector-arrays.
        """

        # The heliostat as a cubic spline ##############################
        hel_spline = self.curve.as_spline(centered=False, rep=False)
        ################################################################

        # The [x,y] point-arrays ####################################################################################
        # Points which define segments of the heliostat surface.
        # Their projected width in the aperture are equal.
        # 'x_values' represent the x-coordinates of central points of segments in the heliostat contour which has
        # a project width on the aperture that is uniform.
        x_coord = array([0.5 * (self.curve.x[i] + self.curve.x[i + 1])
                         for i in range(self.curve.x.shape[0] - 1)])
        # The y-coordinates of the point in the heliostat contour of these 'x_values'.
        y_coord = hel_spline(x_coord)

        seg_pts = zeros(shape=(x_coord.shape[0], 2))
        seg_pts.T[:] = x_coord, y_coord
        ##############################################################################################################

        # The [x,y] normal vectors ###############################
        normals = zeros(shape=(seg_pts.shape[0], 2))
        dy_dx = hel_spline.derivative()

        normals[:] = [nrm(V(arctan(dy_dx(p)))) for p in x_coord]
        normals = R(pi / 2).dot(normals.T).T
        ##########################################################

        # Transforming the [x,y] point and vector-arrays to [x, 0, z] #########
        # seg_points = transform_heliostat(seg_pts).round(3)
        # seg_normals = transform_heliostat(normals).round(6)

        seg_points = seg_pts
        seg_normals = normals
        #######################################################################

        return seg_points, seg_normals

    def as_plane_curve(self):
        return PlaneCurve(curve_pts=self.contour, curve_center=self.center)

    def local_slope(self, weighted=True):

        slope_f = self.curve.as_spline().derivative()
        l_slope = arctan(slope_f(self.seg_pts.T[0])) if weighted else arctan(slope_f(self.contour.T[0]))

        return l_slope


class CylindricalHeliostat:

    def __init__(self,
                 center: array,
                 width: float,
                 theta_d: float = None, aim: array = None,
                 forced_design=False, radius: float = None,
                 nbr_pts=121):

        self.shape = 'cylindrical'
        self.width = abs(width)

        self.center = array([center[0], center[-1]])
        # Ensure an odd number of points
        self.n_pts = nbr_pts + 1 if nbr_pts % 2 == 0 else nbr_pts

        if not forced_design and theta_d is not None and aim is not None:
            self.design_position = theta_d
            self.aim = array([aim[0], aim[-1]])
            self.radius = cylindrical_curvature(center=self.center, aim=self.aim, theta_d=theta_d)
        elif forced_design and radius is not None:
            self.radius = abs(radius)
        else:
            raise ValueError('Invalid arguments. Please, see CylindricalHeliostat class documentation.')

        self.contour = design_cylindrical_mirror(center=self.center, width=self.width, radius=self.radius,
                                                 nbr_pts=self.n_pts)

        # Attribute that holds the Heliostat object as a PlaneCurve object.
        self.curve = self.as_plane_curve()
        self.normals = self.curve.normals2surface()

        self.seg_pts, self.seg_normals = self.segments_of_equal_projected_aperture()

    def segments_of_equal_projected_aperture(self):
        """
        This method calculates [x,0, z] point and normal vectors in the heliostat contour that are the central points
        of segments that has the same projected width in the aperture of the heliostat.
        These points ensure a uniform discretization of the mirror aperture and are essential to compute efficiency
        calculations.

        :return: A tuple of [x,0,z] point-arrays and [x,0,z] vector-arrays.
        """

        # The heliostat as a cubic spline ##############################
        hel_spline = self.curve.as_spline(centered=False, rep=False)
        ################################################################

        # The [x,y] point-arrays ####################################################################################
        # Points which define segments of the heliostat surface.
        # Their projected width in the aperture are equal.
        # 'x_values' represent the x-coordinates of central points of segments in the heliostat contour which has
        # a project width on the aperture that is uniform.
        x_coord = array([0.5 * (self.curve.x[i] + self.curve.x[i + 1])
                         for i in range(self.curve.x.shape[0] - 1)])
        # The y-coordinates of the point in the heliostat contour of these 'x_values'.
        y_coord = hel_spline(x_coord)

        seg_pts = zeros(shape=(x_coord.shape[0], 2))
        seg_pts.T[:] = x_coord, y_coord
        ##############################################################################################################

        # The [x,y] normal vectors ###############################
        normals = zeros(shape=(seg_pts.shape[0], 2))
        dy_dx = hel_spline.derivative()

        normals[:] = [nrm(V(arctan(dy_dx(p)))) for p in x_coord]
        normals = R(pi / 2).dot(normals.T).T
        ##########################################################

        # Transforming the [x,y] point and vector-arrays to [x, 0, z] #########
        # seg_points = transform_heliostat(seg_pts).round(3)
        # seg_normals = transform_heliostat(normals).round(6)
        seg_points = seg_pts
        seg_normals = normals
        #######################################################################

        return seg_points, seg_normals

    def as_plane_curve(self):
        return PlaneCurve(curve_pts=self.contour, curve_center=self.center)

    def local_slope(self, weighted=True):

        slope_f = self.curve.as_spline().derivative()
        l_slope = arctan(slope_f(self.seg_pts.T[0])) if weighted else arctan(slope_f(self.contour.T[0]))

        return l_slope


def design_mirrors(width: float, center: array, aim: array, theta_d: float, nbr_pts: int):

    cyl_mirror = CylindricalHeliostat(width=width,
                                      center=center,
                                      aim=aim,
                                      theta_d=theta_d,
                                      nbr_pts=nbr_pts)

    par_mirror = ParabolicHeliostat(width=width,
                                    center=center,
                                    aim=aim,
                                    theta_d=theta_d,
                                    nbr_pts=nbr_pts)

    return cyl_mirror, par_mirror


def center_from_lambda(lamb: float, h: float):

    lamb_rad = deg2rad(lamb)
    xc = h * tan(lamb_rad)

    return array([xc, 0])


def slope_dev(angle_pair: array, w_ratio, nbr_pts):

    h = 8000
    width = h * w_ratio

    aim = array([0, h])

    lamb, theta_d = angle_pair
    center = center_from_lambda(lamb=lamb, h=h)

    cyl_mirror, par_mirror = design_mirrors(width=width, center=center, aim=aim, theta_d=theta_d, nbr_pts=nbr_pts)

    cyl_slope = cyl_mirror.local_slope(weighted=True)
    par_slope = par_mirror.local_slope(weighted=True)

    flt_slope = zeros(shape=par_slope.shape)

    rms_slope_dev = rmse(predictions=cyl_slope, targets=par_slope) * 1000  # in mrad
    flat_rms_slope_dev = rmse(predictions=flt_slope, targets=par_slope) * 1000  # in mrad

    return rms_slope_dev, flat_rms_slope_dev


def one_ratio_slope_dev(ratio: float, nbr_pts=20):
    angles_list = [[x, y] for x in arange(-70., 72.5, 2.5) for y in arange(-85., 87.5, 2.5)]
    slope_deviations = array([slope_dev(angle_pair=a, w_ratio=ratio, nbr_pts=nbr_pts) for a in angles_list])

    dframe = pd.DataFrame(angles_list, columns=['lambda', 'theta_d'])

    dframe['w_ratio'] = ones(dframe.shape[0]) * ratio
    dframe['n_pts'] = ones(dframe.shape[0]) * nbr_pts

    dframe['rms_slope_dev'] = slope_deviations.T[0]
    dframe['flat_rms_slope_dev'] = slope_deviations.T[1]

    return dframe


def vertical_design_dev(wi: float, we: float, n: int, nbr_pts=401):

    w_range = linspace(start=min(wi, we), stop=max(wi, we), num=n)
    rms_slope_dev = zeros(shape=(w_range.shape[0], 3))

    f = 1
    r = 2
    center = array([0, 0])

    for i, w in enumerate(w_range):

        cyl_mirror = CylindricalHeliostat(width=w, center=center, forced_design=True, radius=r, nbr_pts=nbr_pts)
        par_mirror = ParabolicHeliostat(width=w, center=center, forced_design=True, focal_length=f, nbr_pts=nbr_pts)

        cyl_slope = cyl_mirror.local_slope(weighted=True)
        par_slope = par_mirror.local_slope(weighted=True)
        flt_slope = zeros(shape=par_slope.shape)

        rms_cyl_slop_dev = rmse(predictions=cyl_slope, targets=par_slope) * 1000  # in milliradians
        rms_flt_slop_dev = rmse(predictions=flt_slope, targets=par_slope) * 1000  # in milliradians

        rms_slope_dev[i] = w, rms_cyl_slop_dev, rms_flt_slop_dev

    data = rms_slope_dev

    dframe = pd.DataFrame(data, columns=['w_ratio', 'rms_slope_dev', 'flat_rms_slope_dev'])
    dframe['n_pts'] = ones(dframe.shape[0]) * nbr_pts

    return dframe
