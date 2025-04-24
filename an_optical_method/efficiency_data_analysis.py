from matplotlib import pyplot as plt
from matplotlib.ticker import EngFormatter, PercentFormatter

from solaren.scopy.linear_fresnel import biaxial_annual_eta, annual_eta
from utils import arrays_to_contour, rmse

from an_optical_method.base_models import *
from an_optical_method.efficiency_calculations import files_path as efficiency_files_path

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

# Functions for this program ###########################################################################################


def plot_two_sources_validation(
        lfc: uniform_geometry, sources: (EffectiveSource, EffectiveSource),
        path: Path,
        fig_width: float, fig_height: float, font_size: float,
        fig_path: Path = None, fig_format='pdf', dpi=300.,
        engine='API',
        force_sim=False):

    fig = plt.figure(dpi=dpi, figsize=(fig_width/2.54, fig_height/2.54))

    for i, es in enumerate(sources):

        an_df, rt_df = lfc.validation_analysis(ES=es, path=path, engine=engine,
                                               force_sim=force_sim)

        ax = fig.add_subplot(1, 2, i + 1)

        levels = 13
        x, y, an_z = arrays_to_contour(data=an_df, x_col='theta_t', y_col='theta_l', z_col='eta')
        _, _, rt_z = arrays_to_contour(data=rt_df, x_col='theta_t', y_col='theta_l', z_col='eta')

        cp = ax.contourf(x, y, an_z - rt_z, levels=levels, cmap='viridis')
        ax.contour(x, y, an_z - rt_z, levels=levels, colors='black', linewidths=0.25)

        cb = plt.colorbar(cp)
        cb.set_label(label=r'$\Delta\eta$',
                     fontsize=font_size)

        cb.ax.tick_params(labelsize=font_size-1)

        dev = rmse(predictions=an_df['eta'].values, targets=rt_df['eta'].values).round(4)

        ax.set_xlabel(r'$\theta_T$', fontsize=font_size)
        ax.set_ylabel(r'$\theta_L$', fontsize=font_size)

        plot_title = '(a)' if i == 0 else '(b)'

        ax.set_title(f'{plot_title} LFC {lfc.name[-1]} and {es.name}: ' +
                     r'$\delta_{\eta}$' + f' = {dev}', fontsize=font_size)

        ax.xaxis.set_major_formatter(EngFormatter(unit=r"$^{\circ}$", sep=""))
        ax.yaxis.set_major_formatter(EngFormatter(unit=r"$^{\circ}$", sep=""))
        ax.tick_params(axis='both', which='major', labelsize=font_size-1)

        [spine.set_linewidth(0.35) for spine in ax.spines.values()]

    fig.tight_layout(pad=0.25, w_pad=0.5)
    if fig_path is not None:
        plt.savefig(Path(fig_path,
                         f'validation_{lfc.name}_{sources[0].name}_{sources[1].name}_biaxial.{fig_format}'))
    plt.show()

    return None


def plot_rmse_evolution(
        lfc: uniform_geometry, sources,
        path: Path,
        fig_width: float, fig_height: float,
        font_size=8.,
        fig_path: Path = None, fig_format='pdf', dpi=300.,
        plot_factorized=False,
        force_sim=False,
        engine='API',
        print_source_index=5):

    biaxial_data = zeros(shape=(len(sources), 2))

    tran_data = zeros(shape=(len(sources), 2))
    long_data = zeros(shape=(len(sources), 2))

    for i, es in enumerate(sources):

        an_df, rt_df = lfc.validation_analysis(ES=es, path=path, engine=engine,
                                               force_sim=force_sim)

        biaxial_dev = rmse(predictions=an_df['eta'].values, targets=rt_df['eta'].values).round(4)
        biaxial_data[i] = i + 1, biaxial_dev

        if plot_factorized:
            trans_an = an_df[an_df['theta_l'] == 0.].drop(columns='theta_l').values
            long_an = an_df[an_df['theta_t'] == 0.].drop(columns='theta_t').values

            trans_rt = rt_df[rt_df['theta_l'] == 0.].drop(columns='theta_l').values
            long_rt = rt_df[rt_df['theta_t'] == 0.].drop(columns='theta_t').values

            tran_dev = rmse(predictions=trans_an.T[-1], targets=trans_rt.T[-1]).round(4)
            long_dev = rmse(predictions=long_an.T[-1], targets=long_rt.T[-1]).round(4)

            tran_data[i] = i + 1, tran_dev
            long_data[i] = i + 1, long_dev

    fig = plt.figure(dpi=dpi, figsize=(fig_width/2.54, fig_height/2.54))
    ax = fig.add_subplot()

    ax.plot(*biaxial_data.T, label=r'$\delta_{\eta}$' if plot_factorized else None)
    print(f'Bi-axial RMSE for {effective_sources[print_source_index].name} '
          f'is {biaxial_data[print_source_index].round(4)}')

    if plot_factorized:
        ax.plot(*tran_data.T, label=r'$\delta_{\eta_T}$')
        ax.plot(*long_data.T, label=r'$\delta_{\eta_L}$')
        print(f'Transversal RMSE for {effective_sources[print_source_index].name} '
              f'is {tran_data[print_source_index].round(4)}')
        print(f'Longitudinal RMSE for {effective_sources[print_source_index].name} '
              f'is {long_data[print_source_index].round(4)}')

        ax.legend(fontsize=font_size-1)

    ax.set_xticks([i + 1 for i, es in enumerate(sources)])
    ax.set_xticklabels([es.name for i, es in enumerate(sources)])

    ax.grid(alpha=0.5)
    if plot_factorized:
        ax.set_ylim([0, 1.1*max(biaxial_data.T[-1].max(), tran_data.T[-1].max(), long_data.T[-1].max())])
    else:
        ax.set_ylim([0, 1.1*biaxial_data.T[-1].max()])

    if plot_factorized:
        ax.set_ylabel(r'RMSE', fontsize=font_size)
    else:
        ax.set_ylabel(r'$\delta_{\eta}$', fontsize=font_size)

    ax.tick_params(axis='both', which='major', labelsize=font_size - 1)

    [spine.set_linewidth(0.35) for spine in ax.spines.values()]

    plt.tight_layout(pad=0.25)
    if fig_path is not None:
        plt.savefig(Path(fig_path, f'{lfc.name}_rmse_evolution_by_sources.{fig_format}'))
    plt.show()

    return None


def plot_averaged_eta_per_location(lfc: uniform_geometry, ES: EffectiveSource, locations: tuple,
                                   path: Path,
                                   fig_path: Path,
                                   fig_width: float, fig_height: float, font_size: float,
                                   dpi=300., legend_ax=0,
                                   bar_width=4.0, city_name_pad=-70.,
                                   engine: str = 'API',
                                   from_bottom=False, fig_format='pdf'):

    an_df, rt_df = lfc.validation_analysis(ES=ES, path=path, engine=engine,
                                           force_sim=False)

    latitudes = array([loc.lat for loc in locations])
    names = [loc.name for loc in locations]

    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width/2.54, fig_height/2.54))

    for i in (0, 1):

        ax = fig.add_subplot(1, 2, i + 1)
        ax.set_title(f'(a) NS-mounting' if i == 0 else f'(b) EW-mounting',
                     fontsize=font_size)

        trans_an = an_df[an_df['theta_l'] == 0.].drop(columns='theta_l').values
        long_an = an_df[an_df['theta_t'] == 0.].drop(columns='theta_t').values

        trans_rt = rt_df[rt_df['theta_l'] == 0.].drop(columns='theta_l').values
        long_rt = rt_df[rt_df['theta_t'] == 0.].drop(columns='theta_t').values

        rt_results_data = array([biaxial_annual_eta(biaxial_data=rt_df.values,
                                                    site=loc,
                                                    NS=True if i == 0 else False)
                                 for loc in locations])

        an_results_data = array([biaxial_annual_eta(biaxial_data=an_df.values,
                                                    site=loc,
                                                    NS=True if i == 0 else False)
                                 for loc in locations])

        Delta_eta = 100 * (an_results_data - rt_results_data) / rt_results_data

        rt_results_data_factorized = array([annual_eta(transversal_data=trans_rt,
                                                       longitudinal_data=long_rt,
                                                       site=loc,
                                                       NS=True if i == 0 else False)
                                            for loc in locations])

        an_results_data_factorized = array([annual_eta(transversal_data=trans_an,
                                                       longitudinal_data=long_an,
                                                       site=loc,
                                                       NS=True if i == 0 else False)
                                            for loc in locations])

        Delta_eta_fac = 100 * (an_results_data_factorized - rt_results_data_factorized) / rt_results_data_factorized

        ax.bar(latitudes - bar_width/2, Delta_eta,
               width=bar_width, alpha=0.35, align='center', label='Bi-axial')

        ax.bar(latitudes + bar_width/2, Delta_eta_fac,
               width=bar_width, alpha=0.35, align='center', label='Factorized')

        ax.set_xlabel('Latitude ', fontsize=font_size)
        ax.xaxis.set_major_formatter(EngFormatter(unit=r"$^{\circ}$", sep=""))
        ax.tick_params(labelsize=font_size-1)

        if i == 0:
            ax.set_ylabel(r'$\Delta\bar{\eta}$', fontsize=font_size)
        ax.yaxis.set_major_formatter(PercentFormatter(decimals=1))

        if i == legend_ax:
            ax.legend(fontsize=font_size-1)

        ax2 = ax.twiny()

        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(latitudes)
        ax2.set_xticklabels([name.split()[0] for name in names])

        if from_bottom:
            ax2.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, labelsize=font_size-1)
            ax2.tick_params(axis='x', pad=city_name_pad, labelrotation=90, direction='in', grid_linewidth=0., width=0.)
        else:
            ax2.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True, labelsize=font_size-1)
            ax2.tick_params(axis='x', pad=city_name_pad, labelrotation=90, direction='in', grid_linewidth=0., width=0.)

        ax.grid(alpha=0.2)
        [spine.set_linewidth(0.20) for spine in ax.spines.values()]
        [spine.set_linewidth(0.20) for spine in ax2.spines.values()]

    plt.tight_layout(pad=0.25, w_pad=0.5)
    plt.savefig(Path(fig_path, f'annual_eta_variation_{lfc.name}_{ES.name}.{fig_format}'))
    plt.show()

    return None


########################################################################################################################


figures_path = Path(Path.cwd(), 'efficiency_figures')
figures_path.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':

    geometry = lfc_1
    plot_two_sources_validation(
        lfc=geometry, sources=(ES1, ES5),
        path=efficiency_files_path,
        fig_width=14, fig_height=5, font_size=8,
        fig_path=figures_path)

    plot_rmse_evolution(lfc=geometry, sources=effective_sources,
                        path=efficiency_files_path,
                        plot_factorized=True,
                        fig_width=8, fig_height=5.5, font_size=9,
                        fig_path=figures_path, fig_format='pdf',
                        print_source_index=0)

    plot_averaged_eta_per_location(
        lfc=geometry, ES=ES5, locations=selected_sites,
        path=efficiency_files_path,
        fig_width=12, fig_height=6., font_size=8,
        bar_width=2.,
        city_name_pad=-125,
        fig_path=figures_path, fig_format='pdf')

    plot_averaged_eta_per_location(
        lfc=lfc_2, ES=ES5, locations=selected_sites,
        path=efficiency_files_path,
        fig_width=12, fig_height=6., font_size=8,
        bar_width=2.,
        city_name_pad=-45, legend_ax=1,
        fig_path=figures_path, fig_format='pdf')

    ####################################################################################################################
    ####################################################################################################################
