import seaborn
from matplotlib import pyplot as plt
from numpy import sqrt
from utils import plot_line

from on_secondaries.comparison_Cheng2018.Cheng2018_classes_and_functions import *

########################################################################################################################
# Plots settings #######################################################################################################

seaborn.set(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

########################################################################################################################
########################################################################################################################

# Receiver main data ###################################################################################################

# As from Cheng2018 ######################################################
# reported optimum values -- see Table 3, p. 494, solution 2C.
H_R = 9.44 * 1e3  # in mm

theta_a = 66.29  # degrees
theta_max = 183.69  # degrees
R_ao = 0.1 * 1e3  # in mm. This is absorber tube radius.

# R_OD = 0.11 * 1e3  # in mm. This is the distance between receiver center and CPC vertex.

# other receiver data
# from Table 1, p. 488.
# The distance between the absorber and the glass tube.
# not subject to optimization -- data from the base case: Table 1, p.488.
R_ag = 0.5 * (0.115 - 0.070) * 1e3  # in mm. This is equal to 0.0225 m.

# The outer glass cover radius
R_go = R_ao + R_ag  # in mm


# From this reported optimum values (see Table 3, p. 494. solution 2C),
# it would imply that $\Delta R_{OD}$ would be -0.0125 m, which does not make sense.
# In meters, it follows:

# R_go = R_ao + R_ag => R_go = R_ao + 0.0225. Then,
# R_OD = R_go + DR_OD => R_OD = Rao + 0.0225 + DR_OD
# If R_ao = 0.1 and R_OD = 0.11, then DR_OD is -0.0125, which does not make sense.


# Therefore, here it is assumed that optimization yield the minimum additional gap,
# such that $\Delta R_{OD}$ = 0.01, and the authors wrongly computed $R_{OD} = R_{ao} + \Delta R_{OD}$.
DR_OD = 0.01 * 1e3  # in mm
############################################################################

# Values to be used ########################################################
receiver_height = H_R
cpc_half_acceptance = theta_a
cpc_theta_max = theta_max
secondary_optic_gap = R_go + DR_OD

absorber_tube_radius = R_ao
cover_outer_radius = R_go
cover_inner_radius = cover_outer_radius - 3  # in mm

###########################################################################

########################################################################################################################

# The primary field data ###############################################################################################

# from Table 1, p. 488
number_of_primary_mirrors = 25

# reported optimum values -- see Table 3, p. 494.
primary_mirrors_radius = 24.75 * 1e3  # in mm -- parameter R_m
primary_mirrors_width = 0.25 * 1e3  # in mm -- parameter W_m
primaries_shift_width_ratio = 2  # -- parameter f_{DW}
primary_mirrors_shift = primary_mirrors_width * primaries_shift_width_ratio  # in mm
########################################################################################################################

# Sunshape and optical errors ##########################################################################################

# from text -- first paragraph, Section 4.3.1, p. 493
tracking_error = 0.5 * 1e-3  # in radians
slope_error = 2.5 * 1e-3  # in radians
optical_error = sqrt(4*(tracking_error**2) + 4*(slope_error**2))  # in radians

# from nomenclature table -- Greek symbols, p. 497.
sunshape = 'pillbox'
sun_half_width = 4.65 * 1e-3  # in radians

sun_source = RadialSource(profile=sunshape, size=sun_half_width)

########################################################################################################################

# The Location
cheng_lat, cheng_long = 23.45, 115.90
cheng_loc = SiteData(name='Jiexi County', latitude=cheng_lat, longitude=cheng_lat)

# Optical properties ###################################################################################################

# from Table 1, p. 488
primaries_reflectivity = 0.92
secondary_reflectivity = 0.95
tube_absorptivity = 0.96
glass_cover_transmissivity = 0.95


primaries_reflectance = OpticalProperty.reflector(name='primaries_refl',
                                                  rho=primaries_reflectivity,
                                                  slope_error=0., spec_error=optical_error)

secondary_reflectance = OpticalProperty.secondary(name='secondary_refl',
                                                  rho=secondary_reflectivity,
                                                  slope_error=slope_error, spec_error=0.)

absorber_properties = OpticalProperty.evacuated_tube(alpha=tube_absorptivity,
                                                     tau=glass_cover_transmissivity,
                                                     ref_index=1.52)

optical_properties = OpticalSettings(primaries_property=primaries_reflectance,
                                     secondary_property=secondary_reflectance,
                                     absorber_property=absorber_properties)

########################################################################################################################

contour_points_secondary_optics = 121

# Cheng2018 optic as it is reported ####################################################
cheng_lfc = lfc_optic(name='cheng_optic',
                      nbr_mirrors=number_of_primary_mirrors,
                      mirror_width=primary_mirrors_width,
                      center_distance=primary_mirrors_shift,
                      mirror_radius=primary_mirrors_radius,
                      receiver_height=receiver_height,
                      absorber_radius=absorber_tube_radius,
                      outer_cover_radius=cover_outer_radius,
                      inner_cover_radius=cover_inner_radius)

cheng_lfc.add_oommen_cpc(theta_a=theta_a,
                         theta_max=theta_max,
                         gap_radius=secondary_optic_gap,
                         nbr_contour_points=contour_points_secondary_optics)
#######################################################################################

# An edge-ray CPC for the reported geometry by Chen2018 ###################################
cheng_cpc = lfc_optic(name='cheng_cpc',
                      nbr_mirrors=number_of_primary_mirrors,
                      mirror_width=primary_mirrors_width,
                      center_distance=primary_mirrors_shift,
                      mirror_radius=primary_mirrors_radius,
                      receiver_height=receiver_height,
                      absorber_radius=absorber_tube_radius,
                      outer_cover_radius=cover_outer_radius,
                      inner_cover_radius=cover_inner_radius)

cheng_cpc.add_cpc_secondary(gap_radius=secondary_optic_gap,
                            nbr_contour_points=contour_points_secondary_optics)
###########################################################################################

# An edge-ray CEC for the reported geometry by Chen2018 ###################################
cheng_cec = lfc_optic(name='cheng_cec',
                      nbr_mirrors=number_of_primary_mirrors,
                      mirror_width=primary_mirrors_width,
                      center_distance=primary_mirrors_shift,
                      mirror_radius=primary_mirrors_radius,
                      receiver_height=receiver_height,
                      absorber_radius=absorber_tube_radius,
                      outer_cover_radius=cover_outer_radius,
                      inner_cover_radius=cover_inner_radius)
cheng_cec.add_cec_secondary(gap_radius=secondary_optic_gap,
                            nbr_contour_points=contour_points_secondary_optics)
###########################################################################################


def plot_optic_comparison(
        base_optic: lfc_optic,
        cpc_optic: lfc_optic,
        cec_optic: lfc_optic,
        fig_width: float = 14, fig_height: float = 6, dpi=300,
        font_size: float = 8, line_width: float = 1.,
        ticks_pad=-3, title_pad=3.,
        figure_path: Path = None, fig_format='pdf'):

    lw = abs(line_width)
    t_pad = abs(title_pad)

    t1, t2 = base_optic.tube_edges()

    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width / 2.54, fig_height / 2.54))

    ax1 = fig.add_subplot(1, 2, 1)

    ax1.plot(*plot_line(base_optic.f1 / 1000., t2 / 1000.), color='green', lw=0.5*lw, label='Edge-rays')
    ax1.plot(*plot_line(base_optic.f2 / 1000, t1 / 1000), color='green', lw=0.5*lw)

    [ax1.plot(*hel.T, lw=1.25*lw, label='Primary field' if i == 0 else None, color='black')
     for i, hel in enumerate(base_optic.primary_field.primaries / 1000.)]

    # cheng_lfc.primary_field.plot_primaries(theta_t='horizontal', rec_aim=cheng_lfc.rec_aim,
    #                                        support_size=0.8*cheng_lfc.mirror_width)
    ax1.plot(*base_optic.absorber.absorber_tube.contour.T / 1000,
             color='red', label='Absorber tube', lw=0.5*lw)
    ax1.plot(*base_optic.absorber.outer_tube.contour.T / 1000,
             color='blue', ls='dashed', lw=0.5*lw, label='Glass cover')
    ax1.plot(*base_optic.absorber.inner_tube.contour.T / 1000,
             color='blue', ls='dashed', lw=0.5*lw)

    ax1.legend(fontsize=font_size - 2)

    ax1.set_xlabel('$x$ [m]', fontsize=font_size)
    ax1.set_ylabel('$z$ [m]', fontsize=font_size)
    ax1.tick_params(labelsize=font_size - 1, pad=ticks_pad)

    ax1.set_title('(a) LFC general view', fontsize=font_size, pad=t_pad)

    ax1.grid(lw=0.5*lw)
    [spine.set_linewidth(0.5*lw) for spine in ax1.spines.values()]

    # Plot (b) #########################################################################################################

    ax2 = fig.add_subplot(1, 2, 2)

    ax2.plot(*plot_line(base_optic.f1 / 1000, t2 / 1000), color='green', lw=0.5*lw)
    ax2.plot(*plot_line(base_optic.f2 / 1000, t1 / 1000), color='green', lw=0.5*lw)

    ax2.plot(*base_optic.absorber.absorber_tube.contour.T / 1000, color='red', lw=0.5*lw)
    ax2.plot(*base_optic.absorber.outer_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.5*lw)
    ax2.plot(*base_optic.absorber.inner_tube.contour.T / 1000, color='blue', ls='dashed', lw=0.5*lw)

    ax2.plot(*cpc_optic.secondary.curve_pts.T / 1000, label='CPC', color='black', lw=lw)
    ax2.plot(*cec_optic.secondary.curve_pts.T / 1000, label='CEC', color='brown', lw=lw)
    ax2.plot(*base_optic.secondary.curve_pts.T / 1000, label="Cheng's optic", color='darkorange', lw=lw)

    ax2.legend(ncols=3, loc='upper center', fontsize=font_size - 2)

    ax2.set_xlabel('$x$ [m]', fontsize=font_size)
    ax2.tick_params(labelsize=font_size - 1, pad=ticks_pad)

    ax2.set_title('(b) Receiver view', fontsize=font_size, pad=t_pad)

    ax2.grid(lw=0.5*lw)
    [spine.set_linewidth(0.5*lw) for spine in ax2.spines.values()]

    # Zoom on the receiver view: evacuated tube, secondary optics, and edge-rays
    x0, y0 = base_optic.absorber.center / 1000
    ws = 1.5 * dst(cpc_optic.secondary.curve_pts[0], cpc_optic.secondary.curve_pts[-1]) / 1000
    plt.xlim([x0 - 0.5 * ws, x0 + 0.5 * ws])
    plt.ylim([y0 - 0.75 * ws, y0 + 0.25 * ws])

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'cheng_geometry.{fig_format}'))
    plt.show()

    return None


########################################################################################################################
########################################################################################################################


figures_path = Path(Path.cwd(), 'figures')
figures_path.mkdir(parents=True, exist_ok=True)


if __name__ == '__main__':
    pass

    plot_optic_comparison(
        base_optic=cheng_lfc,
        cpc_optic=cheng_cpc,
        cec_optic=cheng_cec,
        fig_width=13, fig_height=5, font_size=8,
        ticks_pad=-3, title_pad=5.,
        figure_path=figures_path
    )
