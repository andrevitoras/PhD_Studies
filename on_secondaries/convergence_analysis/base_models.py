import seaborn
from matplotlib import pyplot as plt
from utils import plot_line

from on_secondaries.convergence_analysis.classes_and_functions import *

# Configuring plots ################################################################################
seaborn.set_theme(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

####################################################################################################


sun_source = RadialSource(profile='p', size=4.65e-3)

primaries_reflectance = OpticalProperty.reflector(name='primaries_refl',
                                                  rho=0.95, slope_error=2.e-3, spec_error=3.e-3)

secondary_reflectance = OpticalProperty.secondary(name='secondary_refl',
                                                  rho=0.95, slope_error=2.e-3, spec_error=3.e-3)

absorber_properties = OpticalProperty.evacuated_tube(alpha=0.96,
                                                     tau=0.98,
                                                     ref_index=1.52)

optical_properties = OpticalSettings(primaries_property=primaries_reflectance,
                                     secondary_property=secondary_reflectance,
                                     absorber_property=absorber_properties)


mirror_width = 750
center_distance = 1054
nbr_mirrors = 16
receiver_height = 7200
absorber_radius = 35
outer_cover_radius = 62.5
inner_cover_radius = outer_cover_radius - 3.

lfc_1 = lfc_optic(name='novatec_cec',
                  mirror_width=mirror_width, center_distance=center_distance, nbr_mirrors=nbr_mirrors,
                  receiver_height=receiver_height,
                  absorber_radius=absorber_radius,
                  outer_cover_radius=outer_cover_radius,
                  inner_cover_radius=inner_cover_radius,
                  secondary_type='CEC')

lfc_2 = lfc_optic(name='novatec_cpc',
                  mirror_width=mirror_width, center_distance=center_distance, nbr_mirrors=nbr_mirrors,
                  receiver_height=receiver_height,
                  absorber_radius=absorber_radius,
                  outer_cover_radius=outer_cover_radius,
                  inner_cover_radius=inner_cover_radius,
                  secondary_type='CPC')


if __name__ == '__main__':

    t1, t2 = lfc_1.tube_edges()

    fig = plt.figure(dpi=300, figsize=(12, 5))
    ax1 = fig.add_subplot(1, 2, 1)

    lfc_1.primary_field.plot_primaries(theta_t='horizontal', rec_aim=lfc_1.rec_aim, support_size=200)
    ax1.plot(*lfc_1.absorber.absorber_tube.contour.T, color='red')

    ax1.plot(*plot_line(lfc_1.f1, t2), color='green', lw=0.5, label='Edge-rays')
    ax1.plot(*plot_line(lfc_1.f2, t1), color='green', lw=0.5)

    plt.legend(fontsize=10)

    ax2 = fig.add_subplot(1, 2, 2)

    ax2.plot(*lfc_1.secondary.curve_pts.T, label=f'{lfc_1.secondary_type} optic', lw=1)
    ax2.plot(*lfc_2.secondary.curve_pts.T, label=f'{lfc_2.secondary_type} optic', lw=1)

    ax2.plot(*lfc_1.absorber.absorber_tube.contour.T, label='Absorber tube', color='red')
    ax2.plot(*lfc_1.absorber.outer_tube.contour.T, label='Glass cover', color='magenta', ls='dashed', lw=0.5)
    ax2.plot(*lfc_1.absorber.inner_tube.contour.T, color='magenta', ls='dashed', lw=0.5)

    ax2.plot(*plot_line(lfc_1.f1, t2), color='green', lw=0.5)
    ax2.plot(*plot_line(lfc_1.f2, t1), color='green', lw=0.5)

    plt.legend(fontsize=10, ncols=2)

    x0, y0 = lfc_1.absorber.center
    ws = 1.2 * dst(lfc_1.secondary.curve_pts[0], lfc_1.secondary.curve_pts[-1])

    plt.xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    plt.ylim([y0 - 0.75 * ws, y0 + 0.25 * ws])

    plt.show()
