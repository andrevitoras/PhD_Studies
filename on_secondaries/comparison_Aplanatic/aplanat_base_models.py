from copy import deepcopy
from utils import plot_line
from on_secondaries.comparison_Aplanatic.aplanat_functions_and_classes import *

# Configuring plots ####################################################################################################
import seaborn
from matplotlib import pyplot as plt

seaborn.set(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

########################################################################################################################

# Functions for this program ###########################################################################################


def plot_aplanatic_base_model(
        aplanatic_optic: AplanatLfc,
        fig_width: float = 14, fig_height: float = 6, dpi=300,
        font_size: float = 8, line_width: float = 1., marker_size: float = 3,
        title_pad=3.,
        figure_path: Path = None, fig_format='pdf'):

    lw = abs(line_width)
    ms = abs(marker_size)
    t_pad = abs(title_pad)

    number_mirrors = aplanatic_optic.primary_field.nbr_mirrors

    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width/2.54, fig_height/2.54))

    # The first plot -- (a) ############################################################################################
    ax = fig.add_subplot(1, 2, 1)

    [ax.plot(*plot_line(c, a), color='green', lw=0.5)
     for c, a in zip(aplanatic_optic.primary_field.centers / 1000, aplanatic_optic.tracking_points / 1000)]

    ax.plot(*aplanatic_optic.base_aplanat.primary_contour.T / 1000,
            lw=lw, ls='dashed', label='Aplanat primary', color='magenta')
    ax.plot(*aplanatic_optic.base_aplanat.secondary_contour.T / 1000,
            lw=lw, label='Aplanat secondary', color='darkorange')

    ax.plot(*aplanatic_optic.base_aplanat.segment_primary(number_mirrors)[0].T / 1000,
            lw=0, marker='.', ms=ms, color='blue')

    [ax.plot(*hel.T, lw=lw, label='Primary field' if i == 0 else None, color='black')
     for i, hel in enumerate(aplanatic_optic.primary_field.primaries / 1000)]

    ax.plot(*aplanatic_optic.primary_field.centers.T / 1000,
            lw=0, marker='.', ms=0.5*ms, color='red')

    ax.plot(*aplanatic_optic.absorber.contour.T / 1000, lw=lw, label='Absorber tube', color='red')

    ax.tick_params(labelsize=font_size - 1)
    ax.legend(fontsize=font_size - 2, loc='best')

    ax.set_xlabel('$x$ [m]', fontsize=font_size)
    ax.set_ylabel('$z$ [m]', fontsize=font_size)
    ax.set_title(
        f'(a) Base aplanatic optic \n'
        f'$s$ = {aplanatic_optic.base_aplanat.s}, '
        f'$k$ = {aplanatic_optic.base_aplanat.k}, '
        f'NA = {aplanatic_optic.base_aplanat.na}, '
        r'$W_{a}$ = ' + f'{2*aplanatic_optic.base_aplanat.xp_max / 1000} m',
        fontsize=font_size - 1, pad=t_pad)

    # The second plot -- (b) ###########################################################################################

    ax = fig.add_subplot(1, 2, 2)

    [ax.plot(*plot_line(c, a), color='green', lw=0.5*lw)
     for c, a in zip(aplanatic_optic.primary_field.centers / 1000, aplanatic_optic.tracking_points / 1000)]

    ax.plot(*aplanatic_optic.base_aplanat.secondary_contour.T / 1000,
            lw=lw, color='darkorange')
    ax.plot(*aplanatic_optic.absorber.contour.T / 1000, lw=lw, color='red')

    ax.plot(*aplanatic_optic.tracking_points.T / 1000,
            lw=0, marker='.', ms=ms, color='magenta', label='Tracking points')

    Ws = dst(aplanatic_optic.base_aplanat.secondary_contour[0], aplanatic_optic.base_aplanat.secondary_contour[-1])
    Ws = round(Ws / 1000, 3)

    ax.set_xlabel('$x$ [m]', fontsize=font_size)
    ax.set_title('(b) Receiver view \n'
                 f'$r_a$ = {aplanatic_optic.absorber.radius / 1000} m, $W_s$ = {Ws} m',
                 fontsize=font_size - 1, pad=t_pad)

    ax.tick_params(labelsize=font_size - 1)
    ax.legend(fontsize=font_size - 2, loc='best')

    x0, y0 = aplanatic_optic.absorber.center / 1000
    ws = 1.2 * dst(aplanatic_optic.base_aplanat.secondary_contour[0],
                   aplanatic_optic.base_aplanat.secondary_contour[-1]) / 1000

    ax.set_xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    ax.set_ylim([y0 - 0.57 * ws, y0 + 0.43 * ws])

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'aplanatic_linear_fresnel.{fig_format}'))
    plt.show()

    return None


def plot_aplanatic_comparison(
        aplanatic_optic: AplanatLfc,
        aplanatic_cpc: AplanatLfc,
        aplanatic_cec: AplanatLfc,
        fig_width: float = 14, fig_height: float = 6, dpi=300,
        font_size: float = 8, line_width: float = 1., marker_size: float = 3,
        title_pad=3.,
        figure_path: Path = None, fig_format='pdf'):

    lw = abs(line_width)
    ms = abs(marker_size)
    t_pad = abs(title_pad)

    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width / 2.54, fig_height / 2.54))

    t1, t2 = aplanatic_optic.tube_edges()

    # the plot (a) #####################################################################################################
    ax1 = fig.add_subplot(1, 2, 1)

    [ax1.plot(*hel.T, lw=1.5*lw, label='Primary mirrors' if i == 0 else None, color='black')
     for i, hel in enumerate(aplanatic_optic.primary_field.primaries / 1000)]
    ax1.plot(*aplanatic_optic.primary_field.centers.T / 1000, lw=0, marker='.', ms=ms, color='blue')

    ax1.plot(*aplanatic_optic.absorber.contour.T / 1000, lw=lw, label='Absorber tube', color='red')

    ax1.plot(*plot_line(aplanatic_optic.f1 / 1000, t2 / 1000), color='green', lw=0.5*lw, label='Edge-rays')
    ax1.plot(*plot_line(aplanatic_optic.f2 / 1000, t1 / 1000), color='green', lw=0.5*lw)

    ax1.plot(*aplanatic_optic.base_aplanat.secondary_contour.T / 1000, lw=1., color='darkorange')
    ax1.plot(*aplanatic_cpc.secondary.curve_pts.T / 1000, lw=lw, color='black')
    ax1.plot(*aplanatic_cec.secondary.curve_pts.T / 1000, lw=lw, color='brown')

    ax1.set_xlabel('$x$ [m]', fontsize=font_size)
    ax1.set_ylabel('$z$ [m]', fontsize=font_size)
    ax1.set_title('(a) LFC general view', fontsize=font_size, pad=t_pad)

    ax1.legend(fontsize=font_size - 2)
    ax1.tick_params(labelsize=font_size - 1)

    # the plot (b) #####################################################################################################
    ax2 = fig.add_subplot(1, 2, 2)

    ax2.plot(*aplanatic_optic.base_aplanat.secondary_contour.T / 1000, lw=lw, label='Aplanatic', color='darkorange')
    ax2.plot(*aplanatic_cpc.secondary.curve_pts.T / 1000, lw=lw, label='CPC', color='black')
    ax2.plot(*aplanatic_cec.secondary.curve_pts.T / 1000, lw=lw, label='CEC', color='brown')

    ax2.plot(*plot_line(aplanatic_optic.f1 / 1000, t2 / 1000), color='green', lw=lw)
    ax2.plot(*plot_line(aplanatic_optic.f2 / 1000, t1 / 1000), color='green', lw=lw)

    ax2.plot(*aplanatic_optic.absorber.contour.T / 1000, lw=lw, color='red')

    ax2.set_xlabel('$x$ [m]', fontsize=font_size)
    ax2.set_title('(b) Receiver view', fontsize=font_size, pad=t_pad)

    ax2.legend(ncols=3, loc='upper center', fontsize=font_size - 2)
    ax2.tick_params(labelsize=font_size - 1)

    x0, y0 = aplanatic_optic.absorber.center / 1000
    ws = 1.5 * dst(aplanatic_optic.secondary.curve_pts[0], aplanatic_optic.secondary.curve_pts[-1]) / 1000
    ax2.set_xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    ax2.set_ylim([y0 - 0.6 * ws, y0 + 0.4 * ws])

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'aplanatic_geometry.{fig_format}'))
    plt.show()

    return None

########################################################################################################################

# Sunshape and optical errors ##########################################################################################
# optical errors
tracking_error = 0 * 1e-3  # in radians
slope_error = 0 * 1e-3  # in radians
optical_error = 0  # in radians

# sunshape
sunshape = 'gaussian'
sun_size = 9.0 * 1e-3  # in radians
sun_source = RadialSource(profile=sunshape, size=sun_size)
########################################################################################################################

# Optical properties ###################################################################################################

# surface properties
primary_mirrors_reflectance = 0.95
secondary_optic_reflectance = 0.95
absorber_absorbance = 0.95

primaries_reflectance = OpticalProperty.reflector(name='primaries_refl',
                                                  rho=primary_mirrors_reflectance,
                                                  slope_error=0., spec_error=optical_error)

secondary_reflectance = OpticalProperty.secondary(name='secondary_refl',
                                                  rho=secondary_optic_reflectance,
                                                  slope_error=0., spec_error=optical_error)

absorber_properties = OpticalProperty.absorber_tube(name='absorber_abs',
                                                    alpha=absorber_absorbance)

optical_properties = OpticalSettings(primaries_property=primaries_reflectance,
                                     secondary_property=secondary_reflectance,
                                     absorber_property=absorber_properties)

########################################################################################################################


par_s, par_k, par_NA = -2.25, -0.1, 1
aplanat_width = 1.9925 * 1e3
tube_radius = 12

number_of_primary_mirrors = 12

aplanat_lfc = AplanatLfc(name='aplanat_base',
                         s=par_s, k=par_k, NA=par_NA,
                         aplanat_primary_width=aplanat_width,
                         absorber_radius=tube_radius,
                         nbr_primary_mirrors=number_of_primary_mirrors,
                         contour_points=150)

centers_distance = [dst(aplanat_lfc.primary_field.centers[i + 1], aplanat_lfc.primary_field.centers[i])
                    for i in range(len(aplanat_lfc.primary_field.centers) - 1)]

if min(centers_distance) >= aplanat_lfc.primary_field.widths.max():
    print('The shift between mirrors is enough for all primary mirror not collide!')
else:
    print('The shift between mirrors is too short. There are mirrors colliding!')

# CPC and CEC with the gap as the aplanat ###############################
aplanat_lfc_cpc = deepcopy(aplanat_lfc)
aplanat_lfc_cpc.add_cpc_secondary(nbr_contour_points=250)

aplanat_lfc_cec = deepcopy(aplanat_lfc)
aplanat_lfc_cec.add_cec_secondary(nbr_contour_points=250)
#########################################################################

# CPC and CEC with the gap in the proportion of the standard evacuated tube #################
gap_radius = (62.5 / 35) * aplanat_lfc.absorber.radius

aplanat_lfc_cpc_red_gap = deepcopy(aplanat_lfc)
aplanat_lfc_cpc_red_gap.add_cpc_secondary(gap_radius=gap_radius, nbr_contour_points=250)
aplanat_lfc_cpc_red_gap.name = 'aplanat_cpc_red_gap'

aplanat_lfc_cec_red_gap = deepcopy(aplanat_lfc)
aplanat_lfc_cec_red_gap.add_cec_secondary(gap_radius=gap_radius, nbr_contour_points=250)
aplanat_lfc_cec_red_gap.name = 'aplanat_cec_red_gap'
#############################################################################################


figures_path = Path(Path.cwd(), 'figures')
figures_path.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':
    pass

    plot_aplanatic_base_model(
        aplanatic_optic=aplanat_lfc,
        fig_width=14, fig_height=6,
        font_size=8,
        figure_path=figures_path, fig_format='pdf')

    plot_aplanatic_comparison(
        aplanatic_optic=aplanat_lfc,
        aplanatic_cpc=aplanat_lfc_cpc_red_gap,
        aplanatic_cec=aplanat_lfc_cec_red_gap,
        fig_width=14, fig_height=6,
        font_size=8, line_width=0.75, marker_size=3,
        title_pad=5.,
        figure_path=figures_path, fig_format='pdf')
