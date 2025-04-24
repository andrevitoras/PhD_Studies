from decimal import Decimal

import matplotlib.pyplot as plt
import seaborn

from an_optical_method.convergence_analysis.convergence_calculations import *


######################################################
# Plot settings ######################################
seaborn.set_theme(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

########################################################################################################################


# Functions for this program
def plot_analytical_convergence_evaluation(
        lfc: uniform_geometry, ES: EffectiveSource,
        pair_one: (float, float), pair_two: (float, float),
        fig_width: float, fig_height: float, font_size: float,
        fig_name: str = None,
        ticks_pad=-2,
        dpi=300, fig_format='pdf'):

    data_one = convergence_of_points(theta_t=pair_one[0], theta_l=pair_one[1],
                                     base_lfc=lfc, ES=ES, file_path=convergence_files_path, force_sim=False)

    data_two = convergence_of_points(theta_t=pair_two[0], theta_l=pair_two[1],
                                     base_lfc=lfc, ES=ES, file_path=convergence_files_path, force_sim=False)

    fig = plt.figure(dpi=dpi, figsize=(fig_width / 2.54, fig_height / 2.54))
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(*data_one.T)
    ax.set_title(r'(a) \{$\theta_T, \theta_L$\} = \{' +
                 f'{pair_one[0]}' + r'$^{\circ}$' +
                 f',{pair_one[1]}' + r'$^{\circ}$' + '\}', fontsize=font_size)
    ax.set_xlabel(r'$N_p$', fontsize=font_size)
    ax.set_ylabel(r'$\eta^{an}$', fontsize=font_size)

    ax.tick_params(axis='both', which='major', labelsize=font_size-1, pad=ticks_pad)

    ax = fig.add_subplot(1, 2, 2)
    ax.plot(*data_two.T)
    ax.set_title(r'(b) \{$\theta_T, \theta_L$\} = \{' +
                 f'{pair_two[0]}' + r'$^{\circ}$' +
                 f',{pair_two[1]}' + r'$^{\circ}$' + '\}', fontsize=font_size)
    ax.set_xlabel(r'$N_p$', fontsize=font_size)

    ax.tick_params(axis='both', which='major', labelsize=font_size-1, pad=ticks_pad)

    if fig_name is None:
        name = f'an_convergence_{lfc.name}_{ES.name}_{pair_one[0]}_{pair_one[1]}_and_{pair_two[0]}_{pair_two[1]}'
    else:
        name = fig_name

    plt.tight_layout(pad=0.25, w_pad=0.5)
    plt.savefig(Path(convergence_files_path, f'{name}.{fig_format}'))
    plt.show()

    return None


def plot_raytracing_convergence_evaluation(
        lfc: uniform_geometry, ES: EffectiveSource,
        pair_one: (float, float), pair_two: (float, float),
        fig_width: float, fig_height: float, font_size: float,
        fig_name: str = None, line_width=1.,
        ticks_pad=-2,
        dpi=300, fig_format='pdf', engine='API'):

    lw = abs(line_width)

    data_one = convergence_of_ray_hits(theta_t=pair_one[0], theta_l=pair_one[1],
                                       lfc=lfc, base_ES=ES,
                                       file_path=convergence_files_path, force_sim=False, engine=engine)

    data_two = convergence_of_ray_hits(theta_t=pair_two[0], theta_l=pair_two[1],
                                       lfc=lfc, base_ES=ES,
                                       file_path=convergence_files_path, force_sim=False, engine=engine)

    # The first plot
    fig = plt.figure(dpi=dpi, figsize=(fig_width / 2.54, fig_height / 2.54))

    for i, (data, pair) in enumerate(zip([data_one, data_two], [pair_one, pair_two])):

        plot_index = '(a)' if i == 0 else '(b)'
        ax = fig.add_subplot(1, 2, i + 1)

        parts = ax.violinplot([data[k] for k in list(data.keys())],
                              showmeans=True)

        for partname in ('bodies', 'cmins', 'cmaxes', 'cbars'):
            vp = parts[partname]
            if isinstance(vp, list):  # 'bodies' returns a list
                for body in vp:
                    body.set_linewidth(lw)  # Adjust here
            else:
                vp.set_linewidth(lw)

        ax.set_xticks([i + 1 for i in range(len(data.keys()))])
        ax.set_xticklabels([v/(10**6) for v in list(data.keys())])
        ax.tick_params(axis='both', which='major', labelsize=font_size-1, pad=ticks_pad)

        ax.set_title(f'{plot_index} ' + r'\{$\theta_T, \theta_L$\} = \{' +
                     f'{pair[0]}' + r'$^{\circ}$' +
                     f',{pair[1]}' + r'$^{\circ}$' + '\}', fontsize=font_size)

        ax.set_xlabel(r'Ray intersections ($\times 10^6$)', fontsize=font_size)

        if i == 0:
            ax.set_ylabel(r'$\eta^{rt}$', fontsize=font_size)

        ax.grid(alpha=0.35)

    if fig_name is None:
        name = f'{engine}_convergence_{lfc.name}_{ES.name}_{pair_one[0]}_{pair_one[1]}_and_{pair_two[0]}_{pair_two[1]}'
    else:
        name = fig_name

    plt.tight_layout(pad=0.25, w_pad=0.5)
    plt.savefig(Path(convergence_files_path, f'{name}.{fig_format}'))
    plt.show()

    return None


########################################################################################################################


if __name__ == '__main__':

    plot_analytical_convergence_evaluation(
        pair_one=(0, 0), pair_two=(85, 85),
        lfc=lfc_1, ES=ES3,
        fig_width=12, fig_height=5., font_size=8.,
        ticks_pad=-2,
        fig_name='analytical_convergence_analysis')

    plot_raytracing_convergence_evaluation(
        lfc=lfc_1, ES=ES3,
        pair_one=(0, 0), pair_two=(85, 85),
        fig_width=14, fig_height=5.5, font_size=8.,
        line_width=0.75,
        ticks_pad=-2,
        fig_name='raytracing_convergence_analysis')
