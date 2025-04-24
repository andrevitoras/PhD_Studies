from pathlib import Path

import matplotlib.pyplot as plt
import seaborn
from numpy import array, pi

from niopy.geometric_transforms import dst, mid_point
from scopy.linear_fresnel import Absorber
from scopy.linear_fresnel.secondaries import oommen_cpc4tube
from scopy.nio_concentrators import symmetric_cpc2evacuated_tube, oommen_cpc

# Configuring plots #######################################################################
seaborn.set(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

##########################################################################################
"""
This file shows the implementation of Oommen's cpc design [1] as a linear Fresnel secondary optic. For this, it design
and plot the base case of this kind of secondary optic as shown by Qiu et al. [2] and Cheng et al. [3].

References:

[1] Oommen, R., Jayaraman, S., 2001. https://doi.org/10.1016/S0196-8904(00)00113-8.
[2] Qiu et al., 2015. https://doi.org/10.1016/j.apenergy.2015.01.135.
[3] Cheng et al., 2018. A novel optical optimization model for linear Fresnel reflector concentrators.
    https://doi.org/10.1016/j.renene.2018.06.019.

"""


def plot_validation_cpc_optics(
        fig_width: float = 14, fig_height: float = 6,
        font_size=8,
        line_width=1., marker_size=4, grid_alpha=0.75,
        figure_path: Path = None, fig_format='pdf'):

    lw = abs(line_width)
    ms = abs(marker_size)

    # Evacuated absorber data, in meters ################################
    absorber_tube_radius = 0.070 / 2
    outer_cover_radius = 0.115 / 2
    inner_cover_radius = outer_cover_radius - 0.003
    gap_radius = 0.0625
    #####################################################################

    # Qiu [1] CPC optic data ####################################################
    tube_center = array([0, 8.0])  # in meters
    theta_a = 56.0  # in degrees
    theta_max = 3.37 * (180 / pi)  # in degrees

    # Absorber and secondary #####################################################
    qiu_absorber = Absorber.evacuated_tube(name='evacuated_tube',
                                           center=tube_center,
                                           absorber_radius=absorber_tube_radius,
                                           outer_cover_radius=outer_cover_radius,
                                           inner_cover_radius=inner_cover_radius)

    qiu_cpc_optic = oommen_cpc4tube(tube_center=tube_center,
                                    tube_radius=absorber_tube_radius,
                                    gap_radius=gap_radius,
                                    theta_a=theta_a, theta_max=theta_max,
                                    points_per_side=120)

    # geometric data of qiu_cpc_curve
    s1 = qiu_cpc_optic.curve_pts[0]
    s2 = qiu_cpc_optic.curve_pts[-1]
    sm = mid_point(s1, s2)
    aperture_width = dst(s1, s2)
    #############################################################################

    fig = plt.figure(dpi=300,
                     figsize=(fig_width/2.54, fig_height/2.54))

    # Plot (a)
    ax = fig.add_subplot(1, 2, 1)

    # plotting the evacuated absorber tube and optic contour
    ax.plot(*qiu_absorber.absorber_tube.contour.T, label='Absorber', color='red', lw=lw)
    ax.plot(*qiu_cpc_optic.curve_pts.T, label="CPC", color='black', lw=lw)

    # Forcing limits in the axes so that plot present a proper scale for both axis.
    ax.axis('equal')
    ax.grid(alpha=grid_alpha)

    ax.set_xlabel('$x$ [m]',  fontsize=font_size)
    ax.set_ylabel('$z$ [m]',  fontsize=font_size)
    ax.set_title("(a) Qiu's CPC secondary optic",  fontsize=font_size)

    ax.text(x=sm[0], y=sm[1] - 0.025,
            s=r'$W_{s}$ = ' + f'{round(aperture_width, 3)} m\n'
            + r'$h_{s}$ = ' + f'{(tube_center - sm).round(3)[1]} m',
            fontsize=font_size)
    ax.tick_params(labelsize=font_size)
    ax.legend(ncols=2, loc='upper center', fontsize=font_size-1)

    # Plot (b)

    half_acceptance = 60.
    absorber = Absorber.evacuated_tube(name='evacuated_tube',
                                       center=array([0, 0]),
                                       absorber_radius=absorber_tube_radius,
                                       outer_cover_radius=outer_cover_radius,
                                       inner_cover_radius=inner_cover_radius)
    edge_ray_cpc = symmetric_cpc2evacuated_tube(tube_center=array([0, 0]), tube_radius=absorber_tube_radius,
                                                cover_radius=gap_radius, theta_a=half_acceptance,
                                                upwards=True, nbr_pts=10)

    oommen_full_cpc = oommen_cpc(theta_a=half_acceptance, tube_radius=absorber_tube_radius, gap_radius=gap_radius)

    ax = fig.add_subplot(1, 2, 2)
    ax.plot(*absorber.absorber_tube.contour.T, color='red', lw=lw)
    ax.plot(*oommen_full_cpc.T, label="Oommen's design", color='black', lw=lw)
    ax.plot(*edge_ray_cpc[-1].T, label="Chaves' design", color='magenta', lw=0, marker='.', ms=ms)
    ax.plot(*edge_ray_cpc[-2].T, color='magenta', lw=0, marker='.', ms=ms)

    ax.tick_params(labelsize=font_size)
    ax.legend(ncols=1, loc='upper left', fontsize=font_size-1)
    ax.axis('equal')
    ax.grid(alpha=grid_alpha)

    ax.set_xlabel('$x$ [m]', fontsize=font_size)
    ax.set_title("(b) CPCs designs", fontsize=font_size)

    plt.tight_layout(pad=0.25, w_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'nio_designs_validation.{fig_format}'))
    plt.show()

    return None

########################################################################################################################


figures_path = Path(Path.cwd(), 'figures')
figures_path.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':
    pass
    plot_validation_cpc_optics(fig_width=14, fig_height=6,
                               font_size=8, grid_alpha=0.5,
                               figure_path=figures_path, fig_format='pdf')
