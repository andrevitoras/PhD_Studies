from _decimal import Decimal
from numpy import ndarray

from on_secondaries.convergence_analysis.base_models import *


def plot_convergence_analysis(
        simulated_rays: ndarray,
        convergence_results: tuple,
        fig_width: float = 14, fig_height: float = 6, dpi=300,
        font_size: float = 8, line_width: float = 0.5, title_pad=3.,
        figure_path: Path = None, fig_format='pdf'):

    # Preparing data for the plot ##################################################################################
    efficiency_data, uniformity_data, circumferential_uniformity, longitudinal_uniformity = convergence_results
    # data_dic = {'rays_list': simulated_rays.tolist(),
    #
    #             'efficiency': [efficiency_data.mean(axis=1).tolist(), efficiency_data.std(axis=1).tolist()],
    #
    #             'uniformity': [uniformity_data.mean(axis=1).tolist(), uniformity_data.std(axis=1).tolist()],
    #
    #             'circ_uniformity': [circumferential_uniformity.mean(axis=1).tolist(),
    #                                 circumferential_uniformity.std(axis=1).tolist()],
    #
    #             'long_uniformity': [longitudinal_uniformity.mean(axis=1).tolist(),
    #                                 longitudinal_uniformity.std(axis=1).tolist()]
    #             }
    ################################################################################################################

    lw = abs(line_width)
    t_pad = abs(title_pad)

    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width / 2.54, fig_height / 2.54))

    for i, data in enumerate([efficiency_data, circumferential_uniformity]):

        plot_index = '(a)' if i == 0 else '(b)'
        plot_title = 'Optical efficiency' if i == 0 else 'Circumferential uniformity index'

        ax = fig.add_subplot(1, 2, i + 1)

        parts = ax.violinplot([data[k] for k in range(data.shape[0])],
                              showmeans=True)

        for partname in ('bodies', 'cmins', 'cmaxes', 'cbars'):
            vp = parts[partname]
            if isinstance(vp, list):  # 'bodies' returns a list
                for body in vp:
                    body.set_linewidth(lw)  # Adjust here
            else:
                vp.set_linewidth(lw)

        ax.set_xticks([i + 1 for i in range(data.shape[0])])
        ax.set_xticklabels([v/(10**6) for v in simulated_rays])
        ax.tick_params(axis='both', which='major', labelsize=font_size - 1)

        ax.set_title(f'{plot_index} {plot_title}', fontsize=font_size, pad=t_pad)
        ax.set_xlabel(r'Ray intersections ($\times 10^6$)', fontsize=font_size)

        ax.set_ylabel(r'$\eta$' if i == 0 else r'$\delta_q$', fontsize=font_size)
        ax.grid(alpha=0.25)

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'convergence_of_traced_rays.{fig_format}'))
    plt.show()

    return None

########################################################################################################################


simulation_files_path = Path(Path.cwd(), 'simulation_files', 'flux')
simulation_files_path.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':

    sim_per_rays = 50
    rays_to_simulate = array([200e3,
                              1e6, 1.25e6, 1.75e6,
                              2.5e6,
                              3.0e6,
                              4.0e6,
                              5.0e6])

    convergence_data = sensitivity_analysis(
        rays_to_simulate=rays_to_simulate,
        simulations_per_ray=sim_per_rays,
        theta_t=0., dni=1000,
        lfc=lfc_1, optics=optical_properties,
        file_path=simulation_files_path, file_name='optic',
        source=sun_source)

    plot_convergence_analysis(
        simulated_rays=rays_to_simulate,
        convergence_results=convergence_data,
        fig_width=14, fig_height=6,
        font_size=8,
        figure_path=Path.cwd(), fig_format='pdf')


# titles = ['Optical efficiency',
#           'Circumferential uniformity index']
#
# plot_indexes = ['(a)', '(b)']
#
#
# fig = plt.figure(figsize=(12, 5), dpi=300)
# for i, out_data in enumerate([efficiency_data, circumferential_uniformity]):
#
#     ax = fig.add_subplot(1, 2, i + 1)
#
#     plt.plot(rays_to_simulate, out_data.mean(axis=1),
#              lw=1)
#
#     plt.errorbar(rays_to_simulate, out_data.mean(axis=1),
#                  yerr=3*out_data.std(axis=1),
#                  elinewidth=0.5, ecolor='red', capsize=3,
#                  fmt='s', ms=4, color='black')
#
#     plt.title(f'{plot_indexes[i]} {titles[i]}', fontsize=16)
#     plt.xlabel('Number of ray intersections')
#
# plt.tight_layout()
# plt.savefig(f'sensitivity_analysis.svg', bbox_inches='tight', pad_inches=0)
# plt.show()
# ###########################################################
#
# fig = plt.figure(dpi=300, figsize=(12, 5))
# for i, data in enumerate([efficiency_data, circumferential_uniformity]):
#
#     plot_index = '(a)' if i == 0 else '(b)'
#     ax = fig.add_subplot(1, 2, i + 1)
#
#     ax.violinplot([data[k] for k in range(data.shape[0])],
#                   showmeans=True)
#
#     ax.set_xticks([i + 1 for i in range(data.shape[0])])
#     ax.set_xticklabels(["{:.1E}".format(Decimal(v)).replace('+', '') for v in rays_to_simulate])
#     ax.tick_params(axis='both', which='major', labelsize=ticks_font_size)
#
#     ax.set_title(f'{plot_indexes[i]} {titles[i]}', fontsize=16)
#     ax.set_xlabel(r'Ray intersections', fontsize=axes_label_size + 2)
#
#     ax.set_ylabel(r'$\eta$' if i == 0 else r'$\delta_q$')
#     ax.grid(alpha=0.35)
#
# plt.tight_layout()
# plt.savefig(f'convergence_analysis_violin.svg')
# plt.show()
#
