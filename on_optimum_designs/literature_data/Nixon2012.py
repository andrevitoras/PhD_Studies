from niopy.geometric_transforms import dst
from numpy import array, tan, pi, cos, sin, zeros, deg2rad
from scipy.optimize import fsolve
from scopy.linear_fresnel.analysis import shading_loss_factor, tracking_angle


def non_shading_pitch(w: float, xc: float,
                      w_n: float, xc_n: float,
                      tracking_point: array, transversal_incidence: float,
                      in_degrees=True):

    theta_sun = 0.27 * pi / 180.

    theta_t = deg2rad(transversal_incidence) if in_degrees else transversal_incidence
    sm = array([tracking_point[0], 0, tracking_point[-1]])

    hc = array([xc, 0, 0])
    hc_n = array([xc_n, 0, 0])

    tau = abs(tracking_angle(center=hc, rec_aim=sm, theta_t=theta_t, degrees=False))
    tau_n = abs(tracking_angle(center=hc_n, rec_aim=sm, theta_t=theta_t, degrees=False))

    p1 = 0.5 * w * (cos(tau) + sin(tau)*tan(theta_t + theta_sun))
    p2 = 0.5 * w_n * (cos(tau_n) + sin(tau_n)*tan(theta_t + theta_sun))

    return p1 + p2


class nixon_lfc:

    def __init__(self, design_position: float):

        # Units are in meters and degrees
        self.theta_design = abs(design_position)

        self.number_primaries = 28
        self.mirror_width = 0.08

        self.receiver_height = 2.0
        self.nbr_absorber_tubes = 4
        self.absorber_radius = 0.025 / 2

        self.trapezoidal_cavity_aperture = 0.2
        self.trapezoidal_cavity_depth = 0.16
        self.trapezoidal_cavity_back_width = self.nbr_absorber_tubes * 2 * self.absorber_radius

        self.aim_point = array([0, 0, self.receiver_height - self.absorber_radius - self.trapezoidal_cavity_depth])
        # self.aim_point = array([0, 0, self.receiver_height])

        h = self.aim_point[-1]
        theta_s = 0.27 * pi / 180.

        def d(x):

            hc = array([x, 0, 0])
            tau = abs(tracking_angle(center=hc, theta_t=0., rec_aim=self.aim_point))

            wp = 0.5 * self.mirror_width * cos(tau)
            dx = 0.5 * self.mirror_width * sin(tau)
            dx *= tan(theta_s)

            return wp + dx

        d0 = -(0.5 * self.trapezoidal_cavity_aperture + h * tan(theta_s))
        x0 = fsolve(lambda x: d(x[0]) - abs(x[0] - d0), x0=d0 - 0.5*self.mirror_width, full_output=False)[0]

        x_centers = [x0]
        for i in range(self.number_primaries//2 - 1):

            def ds(x):
                s = abs(x - x_centers[i])
                p = non_shading_pitch(w=self.mirror_width, w_n=self.mirror_width,
                                      xc=x_centers[i], xc_n=x,
                                      tracking_point=self.aim_point,
                                      transversal_incidence=self.theta_design, in_degrees=True)

                return s - p

            x_n = fsolve(lambda x: ds(x[0]), x0=x_centers[i] - self.mirror_width, full_output=False)[0]
            x_centers.append(x_n)

        self.x_centers = x_centers

        self.centers = zeros(shape=(len(x_centers), 2))
        self.centers.T[0][:] = x_centers

        self.centers = (-self.centers[::-1]).tolist() + self.centers.tolist()
        self.centers = array(self.centers)

        self.primary_width = dst(self.centers[0], self.centers[-1]) + self.mirror_width

        self.shading_factors = [shading_loss_factor(wi=self.mirror_width, xi=self.x_centers[i],
                                                    wn=self.mirror_width, xn=self.x_centers[i + 1],
                                                    tracking_point=self.aim_point,
                                                    transversal_incidence=self.theta_design, in_degrees=True)
                                for i in range(len(self.x_centers) - 1)]


optic = nixon_lfc(design_position=15)
print(len(optic.centers))
print(round(optic.primary_width, 2))

print(optic.shading_factors)

