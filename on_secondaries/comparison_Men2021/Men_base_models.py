import seaborn
from matplotlib import pyplot as plt
from numpy import sqrt
from utils import plot_line

from on_secondaries.comparison_Men2021.Men_functions_and_classes import *

########################################################################################################################
# Plots settings #######################################################################################################

seaborn.set_theme(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

#######################################################################################################################
########################################################################################################################

# Sunshape and optical errors ##########################################################################################

# optical errors
tracking_error = 0.5 * 1e-3  # in radians
slope_error = 2.5 * 1e-3  # in radians
optical_error = sqrt(4*(tracking_error**2) + 4*(slope_error**2))  # in radians

sunshape = 'pillbox'
sun_half_width = 4.65 * 1e-3  # in radians
sun_source = RadialSource(profile=sunshape, size=sun_half_width)

# Optical properties #####################################

# surface properties
primary_mirrors_reflectance = 0.92
secondary_optic_reflectance = 0.95
absorber_absorbance = 0.96
glass_cover_transmittance = 0.95

primaries_reflectance = OpticalProperty.reflector(name='primaries_refl',
                                                  rho=primary_mirrors_reflectance,
                                                  slope_error=0., spec_error=optical_error)

secondary_reflectance = OpticalProperty.secondary(name='secondary_refl',
                                                  rho=secondary_optic_reflectance,
                                                  slope_error=slope_error, spec_error=0.)

absorber_properties = OpticalProperty.evacuated_tube(alpha=absorber_absorbance,
                                                     tau=glass_cover_transmittance,
                                                     ref_index=1.52)

optical_properties = OpticalSettings(primaries_property=primaries_reflectance,
                                     secondary_property=secondary_reflectance,
                                     absorber_property=absorber_properties)

########################################################################################################################


def plot_men_geometries_comparison(
        index_number: int,
        fig_width: float = 14, fig_height: float = 6, dpi=300,
        font_size: float = 8, line_width: float = 1.,
        ticks_pad=-3, title_pad=3.,
        figure_path: Path = None, fig_format='pdf'):

    base_optic,\
        cpc_optic, cpc_red_gap, \
        cec_optic, cec_red_gap = get_men_optics(index_number=index_number,
                                                secondary_contour_points=contour_points_secondary_optics)

    lw = abs(line_width)
    t_pad = abs(title_pad)

    t1, t2 = base_optic.tube_edges()
    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width / 2.54, fig_height / 2.54))

    ax1 = fig.add_subplot(1, 2, 1)

    l1 = ax1.plot(*plot_line(base_optic.f1 / 1000, t2 / 1000), color='green', lw=0.5*lw, label='Edge-rays')
    ax1.plot(*plot_line(base_optic.f2 / 1000, t1 / 1000), color='green', lw=0.5*lw)

    l2 = ax1.plot(*base_optic.absorber.absorber_tube.contour.T / 1000,
                  color='red', label='Absorber tube', lw=lw)
    l3 = ax1.plot(*base_optic.absorber.outer_tube.contour.T / 1000,
                  color='blue', ls='dashed', lw=0.25*lw, label='Glass cover')
    ax1.plot(*base_optic.absorber.inner_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.25*lw)

    l4 = ax1.plot(*base_optic.secondary.curve_pts.T / 1000,
                  label=f"Men's optic \\#{base_optic.index}", color='darkorange', lw=lw)
    l5 = ax1.plot(*cpc_red_gap.secondary.curve_pts.T / 1000, label='CPC', color='black', lw=0.5*lw)
    l6 = ax1.plot(*cec_red_gap.secondary.curve_pts.T / 1000, label='CEC', color='brown', lw=0.5*lw)

    x0, y0 = base_optic.absorber.center / 1000

    ws = 1.5 * dst(base_optic.secondary.curve_pts[0], base_optic.secondary.curve_pts[-1]) / 1000
    ax1.set_xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    ax1.set_ylim([y0 - 0.65 * ws, y0 + 0.35 * ws])

    lns_low = [l1[0], l2[0], l3[0]]
    lg1 = ax1.legend(handles=lns_low, ncols=3, loc='lower center', fontsize=font_size - 2)
    plt.gca().add_artist(lg1)

    lns_up = [l4[0], l5[0], l6[0]]
    ax1.legend(handles=lns_up, ncols=3, loc='upper center', fontsize=font_size - 2)

    ax1.set_title(f"(a) Men's optic \\#{base_optic.index} with different gap sizes",
                  fontsize=font_size, pad=t_pad)
    ax1.set_xlabel('$x$ [m]', fontsize=font_size)
    ax1.set_ylabel('$z$ [m]', fontsize=font_size)

    ax1.tick_params(labelsize=font_size - 1, pad=ticks_pad)

    # Plot (b) #########################################################################################################
    ax2 = fig.add_subplot(1, 2, 2)

    ax2.plot(*plot_line(base_optic.f1 / 1000, t2 / 1000), color='green', lw=0.5*lw)
    ax2.plot(*plot_line(base_optic.f2 / 1000, t1 / 1000), color='green', lw=0.5*lw)

    ax2.plot(*base_optic.absorber.absorber_tube.contour.T / 1000, color='red', lw=lw)
    ax2.plot(*base_optic.absorber.outer_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.25*lw)
    ax2.plot(*base_optic.absorber.inner_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.25*lw)

    ax2.plot(*base_optic.secondary.curve_pts.T / 1000,
             label=f"Men's optic \\#{base_optic.index}", color='darkorange', lw=lw)
    ax2.plot(*cpc_optic.secondary.curve_pts.T / 1000, label='CPC', color='black', lw=0.5*lw)
    ax2.plot(*cec_optic.secondary.curve_pts.T / 1000, label='CEC', color='brown', lw=0.5*lw)
    # ax2.plot(*cec_optic.secondary.curve_pts[0::4].T / 1000, label='CEC', color='brown', lw=0.5, marker='.', ms=3.)

    ax2.set_yticklabels([])

    ax2.legend(ncols=3, loc='upper center', fontsize=font_size - 2)

    ax2.set_title(f"(b) Men's optic \\#{base_optic.index} with the same gap size",
                  fontsize=font_size, pad=t_pad)
    ax2.set_xlabel('$x$ [m]', fontsize=font_size)

    ax2.tick_params(labelsize=font_size - 1, pad=ticks_pad)

    # General view of the evacuated tube, secondary optics, and edge-rays
    x0, y0 = base_optic.absorber.center / 1000

    ws = 1.5 * dst(base_optic.secondary.curve_pts[0], base_optic.secondary.curve_pts[-1]) / 1000
    ax2.set_xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    ax2.set_ylim([y0 - 0.65 * ws, y0 + 0.35 * ws])

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'men_{base_optic.index}_geometry.{fig_format}'))
    plt.show()

    return None
########################################################################################################################
########################################################################################################################


figures_path = Path(Path.cwd(), 'figures')
figures_path.mkdir(parents=True, exist_ok=True)

contour_points_secondary_optics = 121

########################################################################################################################


if __name__ == '__main__':
    pass

    plot_men_geometries_comparison(
        index_number=15,
        fig_width=14, fig_height=5.5, font_size=7,
        figure_path=figures_path
    )

    # red_gap = True
    # base_optic = men_optic
    #
    # if red_gap:
    #     cpc_optic = men_cpc_red_gap
    #     cec_optic = men_cec_red_gap
    # else:
    #     cpc_optic = men_cpc
    #     cec_optic = men_cec
    #
    # t1, t2 = base_optic.tube_edges()
    # fig = plt.figure(figsize=(12, 5), dpi=300)
    #
    # ax1 = fig.add_subplot(1, 2, 1)
    # plt.title(f"(a) LFC general view: Men's optic \\#{men_optic.index}")
    # plt.xlabel('$x$ [m]')
    # plt.ylabel('$z$ [m]')
    #
    # ax1.plot(*plot_line(base_optic.f1 / 1000., t2 / 1000.), color='green', lw=0.5, label='Edge-rays')
    # ax1.plot(*plot_line(base_optic.f2 / 1000., t1 / 1000.), color='green', lw=0.5)
    #
    # [ax1.plot(*hel.T, lw=1., label='Primary field' if i == 0 else None, color='black')
    #  for i, hel in enumerate(men_optic.primary_field.primaries / 1000.)]
    #
    # ax1.plot(*base_optic.absorber.absorber_tube.contour.T / 1000., color='red', label='Absorber tube')
    # ax1.plot(*base_optic.absorber.outer_tube.contour.T / 1000., color='blue', ls='dashed', lw=0.5, label='Glass cover')
    # ax1.plot(*base_optic.absorber.inner_tube.contour.T / 1000., color='blue', ls='dashed', lw=0.5)
    #
    # plt.legend()
    #
    # ax2 = fig.add_subplot(1, 2, 2)
    # plt.title('(b) Receiver view')
    # plt.xlabel('$x$ [m]')
    #
    # ax2.plot(*plot_line(base_optic.f1 / 1000, t2 / 1000), color='green', lw=0.5)
    # ax2.plot(*plot_line(base_optic.f2 / 1000, t1 / 1000), color='green', lw=0.5)
    #
    # ax2.plot(*base_optic.absorber.absorber_tube.contour.T / 1000, color='red')
    # ax2.plot(*base_optic.absorber.outer_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.5)
    # ax2.plot(*base_optic.absorber.inner_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.5)
    #
    # ax2.plot(*base_optic.secondary.curve_pts.T / 1000,
    #          label=f"Men's optic \\#{men_optic.index}", color='darkorange', lw=1.5)
    # ax2.plot(*cpc_optic.secondary.curve_pts.T / 1000, label='CPC', color='black', lw=0.75)
    # ax2.plot(*cec_optic.secondary.curve_pts.T / 1000, label='CEC', color='brown', lw=0.75)
    # # ax2.plot(*cec_optic.secondary.curve_pts[0::4].T / 1000, label='CEC', color='brown', lw=0.5, marker='.', ms=3.)
    #
    # plt.legend(ncols=3, loc='upper center')
    #
    # # General view of the evacuated tube, secondary optics, and edge-rays
    # x0, y0 = base_optic.absorber.center / 1000
    #
    # ws = 1.5 * dst(base_optic.secondary.curve_pts[0], base_optic.secondary.curve_pts[-1]) / 1000
    # plt.xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    # plt.ylim([y0 - 0.7 * ws, y0 + 0.3 * ws])
    #
    # plt.tight_layout()
    # plt.savefig(Path(figures_path, 'men_nio_geometries.svg'))
    # plt.show()
    #
    # ####################################################################################################################
    #
    # t1, t2 = men_optic.tube_edges()
    # fig = plt.figure(figsize=(12, 5), dpi=300)
    #
    # ax1 = fig.add_subplot(1, 2, 1)
    # plt.title(f"(a) Men's optic \\#{men_optic.index} with different gap sizes")
    # plt.xlabel('$x$ [m]')
    # plt.ylabel('$z$ [m]')
    #
    # l1 = ax1.plot(*plot_line(men_optic.f1 / 1000, t2 / 1000), color='green', lw=0.5, label='Edge-rays')
    # ax1.plot(*plot_line(men_optic.f2 / 1000, t1 / 1000), color='green', lw=0.5)
    #
    # l2 = ax1.plot(*men_optic.absorber.absorber_tube.contour.T / 1000, color='red', label='Absorber tube')
    # l3 = ax1.plot(*men_optic.absorber.outer_tube.contour.T / 1000,
    #               color='blue', ls='dashed', lw=0.5, label='Glass cover')
    # ax1.plot(*men_optic.absorber.inner_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.5)
    #
    # l4 = ax1.plot(*men_optic.secondary.curve_pts.T / 1000,
    #               label=f"Men's optic \\#{men_optic.index}", color='darkorange', lw=1.5)
    # l5 = ax1.plot(*men_cpc_red_gap.secondary.curve_pts.T / 1000, label='CPC', color='black', lw=0.75)
    # l6 = ax1.plot(*men_cec_red_gap.secondary.curve_pts.T / 1000, label='CEC', color='brown', lw=0.75)
    #
    # x0, y0 = men_optic.absorber.center / 1000
    #
    # ws = 1.5 * dst(men_optic.secondary.curve_pts[0], men_optic.secondary.curve_pts[-1]) / 1000
    # plt.xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    # plt.ylim([y0 - 0.65 * ws, y0 + 0.35 * ws])
    #
    # lns_low = [l1[0], l2[0], l3[0]]
    # lg1 = ax1.legend(handles=lns_low, ncols=3, loc='lower center')
    # plt.gca().add_artist(lg1)
    #
    # lns_up = [l4[0], l5[0], l6[0]]
    # ax1.legend(handles=lns_up, ncols=3, loc='upper center')
    #
    # ax2 = fig.add_subplot(1, 2, 2)
    # plt.title(f"(b) Men's optic \\#{men_optic.index} with the same gap size")
    # plt.xlabel('$x$ [m]')
    #
    # ax2.plot(*plot_line(men_optic.f1 / 1000, t2 / 1000), color='green', lw=0.5)
    # ax2.plot(*plot_line(men_optic.f2 / 1000, t1 / 1000), color='green', lw=0.5)
    #
    # ax2.plot(*men_optic.absorber.absorber_tube.contour.T / 1000, color='red')
    # ax2.plot(*men_optic.absorber.outer_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.5)
    # ax2.plot(*men_optic.absorber.inner_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.5)
    #
    # ax2.plot(*men_optic.secondary.curve_pts.T / 1000,
    #          label=f"Men's optic \\#{men_optic.index}", color='darkorange', lw=1.5)
    # ax2.plot(*men_cpc.secondary.curve_pts.T / 1000, label='CPC', color='black', lw=0.75)
    # ax2.plot(*men_cec.secondary.curve_pts.T / 1000, label='CEC', color='brown', lw=0.75)
    # # ax2.plot(*cec_optic.secondary.curve_pts[0::4].T / 1000, label='CEC', color='brown', lw=0.5, marker='.', ms=3.)
    #
    # ax2.set_yticklabels([])
    #
    # plt.legend(ncols=3, loc='upper center')
    #
    # # General view of the evacuated tube, secondary optics, and edge-rays
    # x0, y0 = men_optic.absorber.center / 1000
    #
    # ws = 1.5 * dst(men_optic.secondary.curve_pts[0], men_optic.secondary.curve_pts[-1]) / 1000
    # plt.xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    # plt.ylim([y0 - 0.65 * ws, y0 + 0.35 * ws])
    #
    # plt.tight_layout()
    # plt.savefig(Path(figures_path, 'men_nio_geometries_receiver_view.svg'))
    # plt.show()
    #


