import matplotlib.pyplot as plt
from numpy import linspace

from on_optimum_designs.RQ.base_models import *
from on_optimum_designs.RQ.ea_functions import *

seaborn.set_theme(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'


########################################################################################################################
# Functions for this program ###########################################################################################

reference_labels = {'non-uniform': 'NUN',
                    'variable-radius': 'VAR',
                    'nun-sun-ref': 'NUN-sunRef',
                    'uniform': 'UN',
                    'un-ipa-ref': 'UN-IpaRef'}

reference_colors = {'non-uniform': 'tab:blue',
                    'variable-radius': 'tab:orange',
                    'nun-sun-ref': 'tab:green',
                    'uniform': 'tab:blue',
                    'un-ipa-ref': 'tab:orange'}


def plot_hypervolume_evolution(number_mirrors: int,
                               boundaries: BoundaryConditions,
                               configurations: list,
                               gens: list,
                               labels: list,
                               fig_size=(5, 4),
                               font_size=12):
    fig = plt.figure(dpi=300, figsize=fig_size)
    ax = fig.add_subplot()

    for i, (conf, file, lb) in enumerate(zip(configurations, gens, labels)):
        dic = read_evolved_population_file(file_name=file,
                                           number_of_mirrors=number_mirrors,
                                           configuration=conf,
                                           boundary_conditions=boundaries)

        hv_values = [[j, read_generation_file(file_name=file,
                                              gen=j,
                                              number_of_mirrors=number_mirrors,
                                              configuration=conf,
                                              boundary_conditions=boundaries)['pareto_front']['hv']]

                     for j in range(dic['n_gens'] + 1)]

        hv_values = array(hv_values)

        ax.plot(*hv_values.T, label=lb, lw=0.75)

    ax.legend(fontsize=font_size)
    plt.tight_layout()
    plt.show()

    return None


def plot_comparison(number_mirrors: tuple,
                    boundaries: BoundaryConditions,
                    configurations: list,
                    gen: str,
                    figwidth: float, figheight: float,
                    font_size: float,
                    plot_arrows=True,
                    file_path=None,
                    file_name='configurations_comparison',
                    file_format='png'):

    fig = plt.figure(dpi=300, figsize=(figwidth/2.54, figheight/2.54))
    for i, n in enumerate(number_mirrors):

        min_v = 1000
        max_v = 0

        ax = fig.add_subplot(1, 2, i + 1)
        for j, conf in enumerate(configurations):

            geometric_data = get_pareto_geometries_data(number_of_mirrors=n,
                                                        configuration=conf, boundary_conditions=boundaries,
                                                        gen_name=gen)

            min_v = min(min_v, min(geometric_data['f2']))
            max_v = max(max_v, max(geometric_data['f1']))

            fitness_values = [[f1, f2, ff, sf]
                              for f1, f2, ff, sf in zip(geometric_data['f1'],
                                                        geometric_data['f2'],
                                                        geometric_data['p1'],
                                                        geometric_data['p2'])
                              # if ff < 0.95
                              ]

            fitness_values = array(fitness_values)
            ax.scatter(x=fitness_values.T[0],
                       y=fitness_values.T[1],
                       label=reference_labels[conf], s=1, color=reference_colors[conf])

        if i == 0:
            ax.legend(fontsize=font_size-2)
            ax.set_ylabel('$\Gamma$ [$\mathrm{euros/m^2}$]', fontsize=font_size)
            plot_title = f'(a) $n = $ {n}'
        else:
            plot_title = f'(b) $n = $ {n}'

        ax.set_xlabel('ECF', fontsize=font_size)
        ax.tick_params(axis='both', which='major', labelsize=font_size-1)
        ax.set_title(plot_title, fontsize=font_size)

        ax.axhline(y=min_v, color='red', lw=0.75, ls='dashed')
        ax.axvline(x=max_v, color='red', lw=0.75, ls='dashed')

        if plot_arrows:
            x_lims = ax.get_xlim()
            y_lims = ax.get_ylim()

            ax.annotate(
                r'$\mathrm{ECF}^{\mathrm{max}}$' + f' = {round(max_v, 3)}',  # No text annotation
                xy=(max_v, 0.79 * y_lims[1]),  # End point of the arrow
                xytext=(0.50 * max_v, 0.65 * y_lims[1]),  # Start point of the arrow
                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=4),
                fontsize=font_size-2)

            ax.annotate(
                r'$\Gamma^{\mathrm{min}}$' + f' = {round(min_v, 2)}',  # No text annotation
                xy=(1.2*abs(x_lims[0]), min_v),  # End point of the arrow
                xytext=(1.75*abs(x_lims[0]), min_v*1.5),  # Start point of the arrow
                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=4),
                fontsize=font_size-2)

            # OLD IMPLEMENTATION #######################################################################################
            # if i == 0:
            #     ax.annotate(
            #         f'$f_1$ = {round(max_v, 3)}',  # No text annotation
            #         xy=(max_v, 0.8*y_lims[1]),  # End point of the arrow
            #         xytext=(0.7*max_v, 0.7*y_lims[1]),  # Start point of the arrow
            #         arrowprops=dict(facecolor='black', shrink=0.005, width=2, headwidth=4),
            #         fontsize=font_size-2)
            # else:
            #     ax.annotate(
            #         f'$f_1$ = {round(max_v, 3)}',  # No text annotation
            #         xy=(max_v, 230),  # End point of the arrow
            #         xytext=(max_v - 0.2, 210),  # Start point of the arrow
            #         arrowprops=dict(facecolor='black', shrink=0.005, width=2, headwidth=4),
            #         fontsize=font_size-2)

            # ax.annotate(
            #     f'$f_2$ = {round(min_v, 2)}',  # No text annotation
            #     xy=(0.18, min_v),  # End point of the arrow
            #     xytext=(0.22, min_v + 40),  # Start point of the arrow
            #     arrowprops=dict(facecolor='black', shrink=0.005, width=2, headwidth=4),
            #     fontsize=font_size-2,
            # )
            ############################################################################################################

    plt.tight_layout(pad=0.25, w_pad=0.5)
    if file_path is not None:
        plt.savefig(Path(file_path, f'{file_name}.{file_format}'))
    plt.show()

    return None


def plot_shape_and_filling_factors(number_mirrors: tuple,
                                   boundaries: BoundaryConditions,
                                   configurations: list,
                                   gen: str,
                                   figwidth: float, figheight: float,
                                   font_size: float,
                                   file_path: Path = None,
                                   file_name='shape_and_filling_factor_comparison',
                                   file_format='png'):

    fig = plt.figure(dpi=300,
                     figsize=(figwidth/2.54, figheight/2.54))
    for i, n in enumerate(number_mirrors):

        ax = fig.add_subplot(1, 2, i + 1)
        for j, conf in enumerate(configurations):
            geometric_data = get_pareto_geometries_data(number_of_mirrors=n,
                                                        configuration=conf, boundary_conditions=boundaries,
                                                        gen_name=gen)

            plot_data = [[p1, p2]
                         for p1, p2 in zip(geometric_data['p1'], geometric_data['p2'])]
            plot_data = array(plot_data)

            ax.scatter(*plot_data.T,
                       color=reference_colors[conf],
                       label=reference_labels[conf], s=1)

        if i == 0:
            ax.legend(fontsize=font_size)
            ax.set_ylabel(r'$\pi_2$', fontsize=font_size-1)
            plot_title = f'(a) $n = $ {n}'
        else:
            plot_title = f'(b) $n = $ {n}'

        ax.set_xlabel(r'$\pi_1$', fontsize=font_size)
        ax.tick_params(axis='both', which='major', labelsize=font_size-1)

        ax.set_title(plot_title, fontsize=font_size)

    plt.tight_layout(pad=0.25, w_pad=0.5)
    if file_path is not None:
        plt.savefig(Path(file_path, f'{file_name}.{file_format}'))
    plt.show()

    return None


def plot_parameters_correlation(number_of_mirrors: int,
                                boundaries: BoundaryConditions,
                                configurations: list,
                                gen: str, dni_sum: float,
                                figwidth: float, figheight: float,
                                font_size: float,
                                file_path: Path = None,
                                file_name='parameters_correlation',
                                file_format='png'):

    plot_titles = ['a', 'b', 'c', 'd', 'e', 'f']
    x_labels = [r'$\pi_1$', r'$\pi_2$', r'$\pi_3$']

    x_keys = ['p1', 'p2', 'p3']
    y_keys = ['f1', 'f2']
    keys = [[x, y] for y in y_keys for x in x_keys]

    fig = plt.figure(dpi=300, figsize=(figwidth/2.54, figheight/2.54))
    for i, (x, y) in enumerate(keys):
        ax = fig.add_subplot(3, 3, i + 1)

        for j, conf in enumerate(configurations):
            geometric_data = get_pareto_geometries_data(number_of_mirrors=number_of_mirrors,
                                                        configuration=conf, boundary_conditions=boundaries,
                                                        gen_name=gen)

            ax.scatter(x=geometric_data[x],
                       y=geometric_data[y],
                       label=reference_labels[conf],
                       color=reference_colors[conf],
                       s=0. if x == 'p3' and conf in ['non-uniform', 'nun-sun-ref'] else 1.)

        ax.set_title(f'({plot_titles[i]})', fontsize=font_size)
        ax.tick_params(axis='both', which='major', labelsize=font_size-1)

        if i == 0:
            ax.legend(fontsize=font_size-2)
            ax.set_ylabel(r'ECF', fontsize=font_size)

        if i == 3:
            ax.set_ylabel(r'$\Gamma$ [$\mathrm{euros/m^2}$]', fontsize=font_size)

        # if i < 3:
        #     ax.set_xticklabels([])
        ax.set_xticklabels([])

        if 1 <= i <= 2 or 4 <= i <= 5:
            ax.set_yticklabels([])

        # if i > 2:
        #     ax.set_xlabel(x_labels[i - 3], fontsize=font_size)

        if i == 2 or i == 5:
            x_lim = ax.get_xlim()
            ax.set_xlim([0.95, x_lim[1]])

    plot_titles = ['(g)', '(h)', '(i)']
    for i, x in enumerate(x_keys):

        ax = fig.add_subplot(3, 3, i + 7)
        for j, conf in enumerate(configurations):
            geometric_data = get_pareto_geometries_data(number_of_mirrors=number_of_mirrors,
                                                        configuration=conf, boundary_conditions=boundaries,
                                                        gen_name=gen)
            E = array(geometric_data['f2']) / (dni_sum*array(geometric_data['f1']))
            if i == 0:
                print(f'Minimum value for specific cost of energy is {E.min().round(3)} €/kWh for {conf} configuration')

            ax.scatter(x=geometric_data[x],
                       y=E,
                       label=reference_labels[conf],
                       color=reference_colors[conf],
                       s=0. if x == 'p3' and conf in ['non-uniform', 'nun-sun-ref'] else 1.)

        ax.set_xlabel(x_labels[i], fontsize=font_size)
        ax.set_title(plot_titles[i], fontsize=font_size)
        if i > 0:
            ax.set_yticklabels([])

        ax.tick_params(axis='both', which='major', labelsize=font_size-1)
        if i == 2:
            x_lim = ax.get_xlim()
            ax.set_xlim([0.95, x_lim[1]])

        if i == 0:
            ax.set_ylabel(r'E [euros/kWh]', fontsize=font_size)
            # ax.set_ylabel(r'$\Gamma/\mathrm{ECF} \cdot \sum I_b$ [€/kWh]', fontsize=font_size)

    plt.tight_layout(pad=0.5, h_pad=0.5, w_pad=0.5)
    if file_path is not None:
        plt.savefig(Path(file_path, f'{file_name}.{file_format}'))
    plt.show()

    return None


def plot_variables_frequency_distribution(
        number_of_mirrors: int,
        configurations: list,
        boundaries: BoundaryConditions,
        gen: str,
        figwidth: float, figheight: float,
        num_bins: int = 15,
        font_size: int = 12,
        plot_bounds_lim=False, legend_ax=0,
        file_path=None, file_name='dec_var_freq_distribution', file_format='png'):

    bounds = [boundaries.bounds.height,
              boundaries.bounds.width,
              boundaries.bounds.gap,
              boundaries.bounds.radius]

    fig = plt.figure(dpi=300, figsize=(figwidth/2.54, figheight/2.54))
    axs = [fig.add_subplot(2, 2, i + 1) for i in range(4)]
    x_labels = [r'$H_R$', r'$w$', r'$g$', r'$R$']

    bins_list = []
    [bins_list.append(
        linspace(start=b[0], stop=b[1], num=num_bins) * 1e-3
    )
     for b in bounds]

    for i, conf in enumerate(configurations):
        geometric_data = get_pareto_geometries_data(number_of_mirrors=number_of_mirrors,
                                                    configuration=conf,
                                                    boundary_conditions=boundaries,
                                                    gen_name=gen)

        heights = array(geometric_data['heights']).flatten()

        if conf not in ['non-uniform', 'nun-sun-ref']:
            widths = array(geometric_data['widths']).mean(axis=1)
            gaps = array(geometric_data['gaps']).mean(axis=1)
        else:
            widths = array(geometric_data['widths']).flatten()
            gaps = array(geometric_data['gaps']).flatten()

        if conf in ['uniform', 'un-ipa-ref']:
            radii = array(geometric_data['radii']).mean(axis=1)
        else:
            radii = array(geometric_data['radii']).flatten()

        for j, (ax, values) in enumerate(zip(axs, [heights, widths, gaps, radii])):

            ax.hist(values * 1e-3, alpha=0.25, bins=bins_list[j],
                    color=reference_colors[conf], label=reference_labels[conf])

            # OLD implementation #######################################################################################
            # generate automatically the bins and then use it
            # if i == 0:
            #     n, bins, patches = ax.hist(values * 1e-3, alpha=0.25,
            #                                color=reference_colors[conf], label=reference_labels[conf])
            #     bins_list.append(bins)
            # else:
            #     ax.hist(values * 1e-3, alpha=0.25, bins=bins_list[j],
            #             color=reference_colors[conf], label=reference_labels[conf])
            ############################################################################################################

            ax.grid(alpha=0.15)

            ax.tick_params(axis='both', which='major', labelsize=font_size - 1)
            ax.set_xlabel(f'{x_labels[j]} [m]', fontsize=font_size)

            if j == int(legend_ax):
                ax.legend(fontsize=font_size-1)

            if plot_bounds_lim:
                ax.set_xlim([bounds[j][0] * 1e-3, bounds[j][1] * 1e-3])

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if file_path is not None:
        plt.savefig(Path(file_path, f'n_{number_of_mirrors}_{file_name}.{file_format}'))
    plt.show()

    return None


def plot_var_freq_distribution_per_mirror(
        number_of_mirrors: int,
        configuration: str,
        boundaries: BoundaryConditions,
        gen: str,
        dpi=300, figwidth=12., figheight=12.,
        font_size=12,
        plot_bounds_lim=False, legend_ax=0,
        file_path=None, file_name='var_freq_distribution_per_mirror', file_format='png'):

    # Getting variable bounds to impose axes limits in the plots
    bounds = [boundaries.bounds.height,
              boundaries.bounds.width,
              boundaries.bounds.gap,
              boundaries.bounds.radius]

    # number of mirrors which has decision variables
    nd = number_of_mirrors // 2 if number_of_mirrors % 2 == 0 else 1 + (number_of_mirrors // 2)

    # The figure object and the axes
    fig = plt.figure(dpi=dpi, figsize=(figwidth/2.54, figheight/2.54))
    axs = [fig.add_subplot(nd, 3, i + 1) for i in range(3*nd)]

    geometric_data = get_pareto_geometries_data(number_of_mirrors=number_of_mirrors,
                                                configuration=configuration,
                                                boundary_conditions=boundaries,
                                                gen_name=gen)

    widths = array(geometric_data['widths'])
    gaps = array(geometric_data['gaps'])
    radii = array(geometric_data['radii'])

    bar_alpha = 0.3
    for i in range(nd):

        axs[3*i].hist(widths.T[i] * 1e-3, alpha=bar_alpha,
                      color=reference_colors[configuration], label=reference_labels[configuration])

        axs[3 * i + 1].hist(gaps.T[i] * 1e-3, alpha=bar_alpha,
                            color=reference_colors[configuration], label=reference_labels[configuration])

        axs[3 * i + 2].hist(radii.T[i] * 1e-3, alpha=bar_alpha,
                            color=reference_colors[configuration], label=reference_labels[configuration])

        axs[3 * i].set_ylabel(f'Mirror {i + 1}', fontsize=font_size)

        if i == nd - 1:
            axs[3 * i].set_xlabel(r'$w$ [m]', fontsize=font_size)
            axs[3 * i + 1].set_xlabel(r'$g$ [m]', fontsize=font_size)
            axs[3 * i + 2].set_xlabel(r'$R$ [m]', fontsize=font_size)

    for ax in axs:
        ax.grid(alpha=0.15)
        ax.tick_params(axis='both', which='major', labelsize=font_size - 1)

    plt.tight_layout()
    if file_path is not None:
        plt.savefig(Path(file_path, f'n_{number_of_mirrors}_{file_name}.{file_format}'))
    plt.show()

    return None

########################################################################################################################
