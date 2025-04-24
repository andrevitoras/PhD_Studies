from pathlib import Path

from matplotlib.ticker import EngFormatter
from pandas import concat, read_csv, DataFrame
from scipy.interpolate import interp1d
from scopy.sunlight import RadialSource
from tqdm import tqdm
from utils import closest, arrays_to_contour

from on_the_primaries_shape.functions_and_classes import *

# Plotting configurations ##############################################################################

import matplotlib.pyplot as plt
import seaborn

seaborn.set(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

########################################################################################################
########################################################################################################

# Functions for this program ###########################################################################################


def plot_validation_case(
        w_h_ratio: float,
        lamb: float,
        theta_d: float,
        figure_path: Path,
        fig_width=14., fig_height=8.,
        font_size=10,
        fig_format='pdf'):

    aim = array([0, 8000])
    h = aim[1]
    width = h * w_h_ratio
    center = center_from_lambda(lamb=lamb,
                                h=h)

    cyl_mirror = CylindricalHeliostat(center=center, width=width, theta_d=theta_d, aim=aim,
                                      forced_design=False)
    cyl_mirror_to_plot = CylindricalHeliostat(center=center, width=width, theta_d=theta_d, aim=aim,
                                              forced_design=False, nbr_pts=20)

    par_mirror = ParabolicHeliostat(center=center, width=width, theta_d=theta_d, aim=aim,
                                    forced_design=False)

    # The figure
    fig = plt.figure(dpi=300,
                     figsize=(fig_width/2.54, fig_height/2.54))

    # The first plot shows the mirror contour ###################################
    ax1 = fig.add_subplot(1, 2, 1)

    ax1.set_title(f'(a) Mirrors contour', fontsize=font_size)
    ax1.set_xlabel(r'$x/H$', fontsize=font_size)
    ax1.set_ylabel(r'$z/H$', fontsize=font_size)

    ax1.plot(*par_mirror.contour.T / h, label='Parabolical', lw=1.5)
    ax1.plot(*cyl_mirror_to_plot.contour.T / h, label='Cylindrical', lw=1.5)

    ax1.legend(fontsize=font_size-1)
    ax1.tick_params(axis='both', which='major', labelsize=font_size-1)
    #################################################################################

    # The second plot (b) shows the local slope deviation ###################################
    slop_dev = cyl_mirror.local_slope(weighted=False) - par_mirror.local_slope(weighted=False)
    rms_dev = 1000 * rmse(predictions=cyl_mirror.local_slope(weighted=True),
                          targets=par_mirror.local_slope(weighted=True))

    ax2 = fig.add_subplot(1, 2, 2)

    ax2.plot(cyl_mirror.contour.T[0] / h,
             slop_dev, color='black', lw=1)

    ax2.plot(cyl_mirror.contour.T[0] / h, slop_dev,
             lw=0,
             label=r"$\delta_{\xi}$ = " + f'{round(rms_dev, 1)} mrad')

    ax2.legend(fontsize=font_size, loc='upper center')
    ax2.set_title(r'(b) Local slope deviation',
                  fontsize=font_size)
    ax2.set_xlabel(r'$x/H$', fontsize=font_size)
    ax2.set_ylabel(r'$\Delta \xi$ [rad]', fontsize=font_size)

    ax2.tick_params(axis='both', labelsize=font_size-1)
    #################################################################################

    # Saving figure
    fig.tight_layout(pad=0.25, w_pad=1.)
    plt.savefig(Path(figure_path, f'validation_case.{fig_format}'))
    plt.show()

    return None


def plot_Cheng_Qiu_cases(fig_width=8., fig_height=6.,
                         font_size=9.,
                         figure_path: Path = None, fig_format='pdf'):

    center = array([0, 0, 0])
    nbr_pts = 401

    Qiu_c = CylindricalHeliostat(width=600,
                                 radius=21.7 * 10 ** 3,
                                 forced_design=True,
                                 center=center,
                                 nbr_pts=nbr_pts)

    Qiu_p = ParabolicHeliostat(width=600,
                               focal_length=10.6 * 10 ** 3,
                               center=center,
                               nbr_pts=nbr_pts,
                               forced_design=True)

    Cheng_c = CylindricalHeliostat(width=250,
                                   radius=24.75 * 10 ** 3,
                                   forced_design=True,
                                   center=center,
                                   nbr_pts=nbr_pts)

    Cheng_p = ParabolicHeliostat(width=250,
                                 focal_length=12.33 * 10 ** 3,
                                 center=center,
                                 nbr_pts=nbr_pts,
                                 forced_design=True)

    Qiu_slope_dev = Qiu_c.local_slope(weighted=False) - Qiu_p.local_slope(weighted=False)
    Cheng_slope_dev = Cheng_c.local_slope(weighted=False) - Cheng_p.local_slope(weighted=False)

    Qiu_slope_dev_RMS = 1000 * rmse(predictions=Qiu_c.local_slope(), targets=Qiu_p.local_slope())
    Cheng_slope_dev_RMS = 1000 * rmse(predictions=Cheng_c.local_slope(), targets=Cheng_p.local_slope())

    qiu_rms_value = Qiu_slope_dev_RMS.round(2)
    cheng_rms_value = Cheng_slope_dev_RMS.round(2)

    fig = plt.figure(dpi=300,
                     figsize=(fig_width/2.54, fig_height/2.54))
    ax = fig.add_subplot()
    ax.plot(Qiu_c.contour.T[0] / 1000, Qiu_slope_dev,
            label='Qiu et al. [6] (' + r'$\delta_{\xi} = $' + f' {qiu_rms_value} mrad)')

    ax.plot(Cheng_c.contour.T[0] / 1000, Cheng_slope_dev,
            label='Cheng et al. [8] (' + r'$\delta_{\xi} = $' + f' {cheng_rms_value} mrad)')

    ax.set_xlabel(r'$x$ [m]', fontsize=font_size)
    ax.set_ylabel(r'$\Delta \xi$ [rad]',fontsize=font_size)
    ax.tick_params(labelsize=font_size-1)
    ax.legend(fontsize=font_size-1)

    fig.tight_layout(pad=0.25, w_pad=1.)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'Qiu_and_Cheng_deviations.{fig_format}'))
    plt.show()

    return None


def plot_worst_case_design_domain(dev_database: DataFrame,
                                  fig_width=8, fig_height=8,
                                  font_size=9,
                                  figure_path: Path = None, file_format='pdf'):

    df = dev_database[dev_database['w_ratio'] <= 0.5]
    w_ratios = array(list(set(df['w_ratio'])))
    w_ratios.sort()

    max_dev = zeros(w_ratios.shape[0])
    max_dev_flat = zeros(w_ratios.shape[0])

    for i, w in enumerate(w_ratios):
        data = df[df['w_ratio'] == w]
        max_dev[i] = data['rms_slope_dev'].max()
        max_dev_flat[i] = data['flat_rms_slope_dev'].max()

    fig = plt.figure(dpi=300,
                     figsize=(fig_width/2.54, fig_height/2.54))
    ax = fig.add_subplot()

    ax.plot(w_ratios, max_dev)

    # ax.plot(w_ratios, max_dev, label='Cylindrical mirrors')
    # ax.plot(w_ratios, max_dev_flat, label='Flat mirrors')
    # ax.legend(fontsize=font_size-1)

    ax.set_xlabel(r'$w/H$', fontsize=font_size)
    ax.set_ylabel(r'$\delta_{\xi,max}$ [mrad]', fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size-1)

    fig.tight_layout(pad=0.25)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'max_dev_whole_design_domain.{file_format}'))
    plt.show()

    return None


def plot_worst_case_specific_reference(dev_database: DataFrame,
                                       fig_width=9.5, fig_height=7,
                                       font_size=9,
                                       figure_path: Path = None, fig_format='pdf'):

    database = dev_database
    specific_ref_function = interp1d(database['w_ratio'], database['rms_slope_dev'])

    Abbas2012_ratio = 0.6 / 8
    delta_xi_Abbas2012 = specific_ref_function(Abbas2012_ratio).round(4)

    Boito2017_ratio = 0.2
    delta_xi_Boito2017 = specific_ref_function(Boito2017_ratio).round(4)

    Balaji2016_ratio = round(1.07 / 8, 3)
    delta_xi_Balaji2016 = specific_ref_function(Balaji2016_ratio).round(4)

    fig = plt.figure(dpi=300,
                     figsize=(fig_width/2.54, fig_height/2.54))
    ax = fig.add_subplot()

    ax.plot(database['w_ratio'], database['rms_slope_dev'])

    ax.plot([Abbas2012_ratio], delta_xi_Abbas2012, label='Abbas et al. [4]',
            linewidth=0, marker='.', ms=8, color='red')

    ax.plot(Balaji2016_ratio, delta_xi_Balaji2016, label='Balaji et al [5]',
            linewidth=0, marker='.', ms=8)

    ax.plot(Boito2017_ratio, delta_xi_Boito2017, label='Boito and Grena [7]',
            linewidth=0, marker='.', ms=8, color='green')

    ax.annotate(f'({Abbas2012_ratio}, {delta_xi_Abbas2012})', fontsize=font_size - 1,
                xy=(Abbas2012_ratio, delta_xi_Abbas2012), xytext=(Abbas2012_ratio, 0.16),
                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=5))

    ax.annotate(f'({Balaji2016_ratio}, {delta_xi_Balaji2016.round(3)})', fontsize=font_size - 1,
                xy=(Balaji2016_ratio, delta_xi_Balaji2016), xytext=(Balaji2016_ratio, 0.12),
                arrowprops=dict(facecolor='black', shrink=0.075, width=2, headwidth=5))

    ax.annotate(f'({Boito2017_ratio}, {delta_xi_Boito2017.round(3)})', fontsize=font_size - 1,
                xy=(Boito2017_ratio, delta_xi_Boito2017), xytext=(Boito2017_ratio, 0.09),
                arrowprops=dict(facecolor='black', shrink=0.09, width=2.0, headwidth=5))

    ax.legend(fontsize=font_size - 1)

    ax.set_xlabel(r'$w/H$', fontsize=font_size)
    ax.set_ylabel(r'$\delta_{\xi,max}$ [mrad]', fontsize=font_size)

    ax.tick_params(labelsize=font_size-1)

    fig.tight_layout(pad=0.25)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'slope_dev_spec_ref_design.{fig_format}'))
    plt.show()

    return None


def plot_design_tools(dev_database: DataFrame,
                      sp_ref_database: DataFrame,
                      fig_width=10, fig_height=7,
                      font_size=9,
                      figure_path: Path = None, fig_format='pdf'):

    # Worst case scenario ####################################################################################
    df = dev_database[dev_database['w_ratio'] <= 0.5]
    w_ratios = array(list(set(df['w_ratio'])))
    w_ratios.sort()

    max_dev = zeros(w_ratios.shape[0])
    max_dev_flat = zeros(w_ratios.shape[0])

    for i, w in enumerate(w_ratios):
        data = df[df['w_ratio'] == w]
        max_dev[i] = data['rms_slope_dev'].max()
        max_dev_flat[i] = data['flat_rms_slope_dev'].max()

    max_slope_dev_function = interp1d(w_ratios, max_dev, kind='cubic')
    ########################################################################################################

    specific_ref_function = interp1d(sp_ref_database['w_ratio'], sp_ref_database['rms_slope_dev'])

    sun_csr005 = RadialSource(profile='b', size=0.05)
    opt_error = RadialSource(profile='g', size=5.0)
    source_rms_width = (sun_csr005.rms_width**2 + opt_error.rms_width**2)**0.5

    database = dev_database
    df = database[(database['theta_d'] == 0)]

    w_ratios = array(list(set(df['w_ratio'])))
    w_ratios.sort()

    zen_max_dev = zeros(w_ratios.shape[0])
    for i, w in enumerate(w_ratios):
        data = df[df['w_ratio'] == w]
        zen_max_dev[i] = data['rms_slope_dev'].max()
    zen_dev_function = interp1d(w_ratios, zen_max_dev)

    ratios_to_plot = linspace(start=w_ratios[0], stop=0.2, num=50)

    fig = plt.figure(dpi=300,
                     figsize=(fig_width/2.54, fig_height/2.54))
    ax = fig.add_subplot(1, 2, 1)

    ax.plot(ratios_to_plot, max_slope_dev_function(ratios_to_plot), lw=1, label='All design domain')
    ax.plot(ratios_to_plot, zen_dev_function(ratios_to_plot), lw=1, label='Zenithal reference')
    ax.plot(ratios_to_plot, specific_ref_function(ratios_to_plot), lw=1, label='Specific reference')

    ax.set_xlabel(r'$w/H$', fontsize=font_size)
    ax.set_ylabel(r'$\delta_{\xi,max}$ [mrad]', fontsize=font_size)
    ax.legend(fontsize=font_size - 1)
    ax.tick_params(labelsize=font_size - 1)
    ax.set_title(f'(a)', fontsize=font_size)

    ax = fig.add_subplot(1, 2, 2)

    ax.plot(ratios_to_plot,
            100 * (-1 + (((source_rms_width**2 + 4*max_slope_dev_function(ratios_to_plot)**2)**0.5) / source_rms_width)),
            lw=1, label='All design domain')
    ax.plot(ratios_to_plot,
            100 * (-1 + (((source_rms_width**2 + 4*zen_dev_function(ratios_to_plot)**2)**0.5) / source_rms_width)),
            lw=1, label='Zenithal reference')
    ax.plot(ratios_to_plot,
            100 * (-1 + (((source_rms_width**2 + 4*specific_ref_function(ratios_to_plot)**2)**0.5) / source_rms_width)),
            lw=1, label='Specific reference')

    ax.set_xlabel(r'$w/H$', fontsize=font_size)
    ax.set_ylabel(r'$\Delta\delta_{es}$', fontsize=font_size)
    ax.yaxis.set_major_formatter(EngFormatter(unit=r"\%", sep="", places=1))
    ax.set_yticks([0, 1., 2., 3., 4.])
    ax.tick_params(labelsize=font_size - 1)
    ax.set_title(f'(b)', fontsize=font_size)

    fig.tight_layout(pad=0.25, w_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'design_tools.{fig_format}'))
    plt.show()

    print(f'Maximum deviation for w/H = 0.2 is {max_slope_dev_function(0.2).round(3)}')

    print(f'A Buie CSR 5% has an RMS width of {round(sun_csr005.rms_width, 3)} mrad.')
    print(f'An Gaussian optical error with 5.0 mrad '
          f'standard deviation has an RMS width of {round(opt_error.rms_width, 3)} mrad.')

    print(f'The effective source has an RMS width of {round(source_rms_width, 3)} mrad')

    return None


def plot_design_domain_cases(dev_database: DataFrame,
                             ratio_cases=(0.1, 0.3),
                             fig_width=14., fig_height=7.0,
                             font_size=8.,
                             figure_path: Path = None, fig_format='pdf'):

    seaborn.set(style='white')
    plt.rc('text', usetex=True)
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'NewComputerModern10'

    df = dev_database
    w_ratios = array(list(set(df['w_ratio'])))
    w_ratios.sort()

    ratio_values = [closest(list(w_ratios), ratio) for ratio in ratio_cases]

    fig = plt.figure(dpi=300,
                     figsize=(fig_width/2.54, fig_height/2.54))

    for i, r in enumerate(ratio_values):
        ax = fig.add_subplot(1, 2, i + 1)
        ax.tick_params(axis='both', which='major', labelsize=font_size - 1)

        plot_index = '(a)' if i == 0 else '(b)'
        color_map = 'cividis' if i == 1 else 'viridis'

        x, y, z = arrays_to_contour(data=df[df['w_ratio'] == r],
                                    x_col='lambda',
                                    y_col='theta_d',
                                    z_col='rms_slope_dev')

        cp = ax.contourf(x, y, z, levels=15, cmap=color_map)
        cb = plt.colorbar(cp)
        cb.set_label(label=r'$\delta_{\xi}$ [mrad]', fontsize=font_size)
        cb.ax.tick_params(labelsize=font_size - 1)

        ax.contour(x, y, z, levels=15, colors='black', linewidths=0.5)
        ax.axline(xy1=(-70, -70), xy2=(70, 70), color='red', lw=1)

        ax.set_xlabel(r'$\lambda$', fontsize=font_size)
        ax.set_ylabel(r'$\theta_d$', fontsize=font_size)
        ax.set_title(f'{plot_index} $w/H=$ {round(r, 1)} \n(' + r'$\delta_{\xi ,max} = $'
                     + f' {z.max().round(2)} mrad)',
                     fontsize=font_size)

        ax.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
        ax.yaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))

    fig.tight_layout(pad=0.25, w_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'slope_dev_two_cases.{fig_format}'))
    plt.show()

    return None

###########################################################################################


if __name__ == '__main__':

    # Calculations #####################################################################################################
    ratios = linspace(start=40 / 8000, stop=4000 / 8000, num=50)

    # All design positions
    file_name = 'slope_dev_database'
    file_path = Path.cwd()
    full_database = Path(file_path, f"{file_name}.csv")

    if full_database.is_file():
        full_df = read_csv(full_database)
    else:
        frames = [one_ratio_slope_dev(ratio=r, nbr_pts=401) for r in tqdm(ratios)]
        full_df = concat(frames, ignore_index=True)
        full_df.to_csv(full_database)

    # Only specific reference
    file_name = 'specific_ref_dev_database'
    file_path = Path.cwd()
    specific_ref_database = Path(file_path, f"{file_name}.csv")

    if specific_ref_database.is_file():
        sr_df = read_csv(specific_ref_database)
    else:
        sr_df = vertical_design_dev(wi=ratios[0], we=ratios[-1], n=ratios.shape[0], nbr_pts=401)
        sr_df.to_csv(specific_ref_database)

    ####################################################################################################################

    # Plots #####################################################################################################
    figures_path = Path(Path.cwd(), 'Figures')
    figures_path.mkdir(parents=True, exist_ok=True)

    plot_validation_case(w_h_ratio=0.3, lamb=20, theta_d=-50,
                         figure_path=figures_path,
                         fig_width=13, fig_height=6, font_size=8)

    plot_worst_case_design_domain(dev_database=full_df,
                                  fig_width=8, fig_height=6,
                                  font_size=9,
                                  figure_path=figures_path)

    plot_Cheng_Qiu_cases(fig_width=8, fig_height=6,
                         font_size=9,
                         figure_path=figures_path)

    plot_worst_case_specific_reference(dev_database=sr_df,
                                       fig_width=10, fig_height=7,
                                       font_size=9,
                                       figure_path=figures_path)

    plot_design_tools(dev_database=full_df, sp_ref_database=sr_df,
                      fig_width=12, fig_height=6,
                      font_size=9,
                      figure_path=figures_path)

    plot_design_domain_cases(dev_database=full_df,
                             fig_width=14, fig_height=6,
                             font_size=9,
                             figure_path=figures_path)
