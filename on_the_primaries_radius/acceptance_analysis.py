from pathlib import Path

from matplotlib.ticker import EngFormatter, StrMethodFormatter

from on_the_primaries_radius.classes_and_functions import acceptance4designs, avg_acc4designs
import seaborn
from matplotlib import pyplot as plt
from numpy import array, arange

from solaren.scopy.sunlight import RadialSource, SiteData
from on_the_primaries_radius.base_models import evora, aswan, lfc_1, lfc_2, ES1, ES2, lfr_geometry

######################################################
######################################################


######################################################
# Plot settings ######################################

########################################################################################################################
# Functions for this program  ##########################################################################################

# OLD FUNCTIONS ########################################################################################################
# def plot_non_uniform_acceptances(geometry: lfr_geometry, location: SiteData, sources: tuple,
#                                  figure_path: Path = None,
#                                  line_width=1.25,
#                                  fig_size=(12, 8), dpi=450,
#                                  legend_font_size=14, ticks_font_size=14, axes_label_size=18):
#
#     assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'
#
#     seaborn.set(style='whitegrid')
#     seaborn.set_context('notebook')
#
#     plt.rc('text', usetex=True)
#     plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
#     plt.rcParams['font.family'] = 'serif'
#     plt.rcParams['font.serif'] = 'Latin Modern Roman'
#     plt.rcParams['xtick.labelsize'] = ticks_font_size
#     plt.rcParams['ytick.labelsize'] = ticks_font_size
#
#     lw = abs(line_width)
#
#     fig = plt.figure(figsize=fig_size, dpi=dpi)
#     fig.suptitle(
#         f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$' + f' {round(location.lat, 2)}' + r'$^{\circ})$',
#         fontsize=axes_label_size)
#
#     for i, es in enumerate(sources):
#
#         data = read_acceptance(site=location, geom=geometry, source=es)
#         ns, ew = data
#
#         ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
#         ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)
#
#         for k, dic in enumerate([ns, ew]):
#
#             ax = ax1 if k == 0 else ax2
#
#             ax.plot(*dic['sun'].T, color='green', label='Sun reference', lw=lw)
#             ax.plot(*dic['spr'].T, color='red', label='Specific reference', lw=lw)
#             ax.plot(*dic['boito'].T, color='blue', label='BG design', lw=lw)
#             ax.plot(*dic['nun'].T, color='magenta', label='NUN-OR', lw=lw)
#
#             if i == 0:
#                 plot_index = '(a)' if k == 0 else '(b)'
#
#             else:
#                 plot_index = '(c)' if k == 0 else '(d)'
#
#             ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
#                          fontsize=legend_font_size + 2)
#
#             if k == 0:
#                 ax.set_ylabel(r"$\bar{\eta}$", fontsize=axes_label_size)
#
#             if i == 0:
#                 ax.xaxis.set_ticklabels([])
#             else:
#                 ax.set_xlabel(r"$\theta_{T}$", fontsize=axes_label_size)
#                 ax.xaxis.set_major_formatter(EngFormatter(unit=r"$^{\circ}$", sep=""))
#
#         if i == 0:
#             ax1.legend(fontsize=legend_font_size)
#
#     plt.tight_layout()
#     if figure_path is not None:
#         plt.savefig(Path(figure_path, f'{geometry.name}_{location.name}_non_uniform_acceptances.svg'))
#     plt.show()
#
#     return None
#
# def plot_non_uniform_avg_acceptances(geometry: lfr_geometry, location: SiteData, sources: tuple,
#                                      figure_path: Path = None,
#                                      fig_size=(10, 6), dpi=300,
#                                      ticks_font_size=14, axes_label_size=14):
#
#     seaborn.set(style='darkgrid')
#     seaborn.set_context('notebook')
#
#     plt.rc('text', usetex=True)
#     plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
#     plt.rcParams['font.family'] = 'serif'
#     plt.rcParams['font.serif'] = 'Latin Modern Roman'
#     plt.rcParams['xtick.labelsize'] = ticks_font_size
#     plt.rcParams['ytick.labelsize'] = ticks_font_size
#
#     dic_keys = ['sun', 'spr', 'boito', 'nun']
#     bar_labels = ['Sun\nreference', 'Specific\nreference', 'BG design', 'NUN-OR']
#
#     fig = plt.figure(figsize=fig_size, dpi=dpi)
#     fig.suptitle(
#         f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$' + f' {round(location.lat, 2)}' + r'$^{\circ})$',
#         fontsize=axes_label_size)
#
#     for i, es in enumerate(sources):
#
#         avg_acceptance_data = avg_acc4designs(geom=geometry, site=location, source=es)
#         ns, ew = avg_acceptance_data
#
#         ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
#         ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)
#
#         for k, dic in enumerate([ns, ew]):
#
#             ax = ax1 if k == 0 else ax2
#
#             heights = [dic[k] for k in dic_keys]
#             heights = array(heights).round(2)
#             bar_pos = arange(len(bar_labels))
#
#             ax.bar(bar_pos, heights, align='center')
#             add_labels(bar_pos, heights,
#                        ax=ax,
#                        font_size=ticks_font_size)
#
#             # Create names on the x-axis
#             ax.set_ylim([0, 1.12*heights.max()])
#             ax.yaxis.set_ticklabels([])
#
#             if k == 0:
#                 ax.set_ylabel(r'$\bar{\theta}_a$', fontsize=axes_label_size)
#
#             if i == 1:
#                 ax.set_xticks(bar_pos, bar_labels, rotation=0, ha='center', fontsize=axes_label_size)
#             else:
#                 ax.set_xticks([])
#
#             if i == 0:
#                 plot_index = '(a)' if k == 0 else '(b)'
#
#             else:
#                 plot_index = '(c)' if k == 0 else '(d)'
#
#             ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
#                          fontsize=axes_label_size)
#
#     plt.tight_layout()
#     if figure_path is not None:
#         plt.savefig(Path(figure_path, f'{geometry.name}_{location.name}_non_uniform_avg_acceptances.svg'))
#     plt.show()
#
#     return None
#
# def plot_uniform_avg_acceptances(geometry: lfr_geometry, location: SiteData, sources: tuple,
#                                  figure_path: Path = None,
#                                  fig_size=(10, 6), dpi=300,
#                                  ticks_font_size=14, axes_label_size=14):
#
#     seaborn.set(style='darkgrid')
#     seaborn.set_context('notebook')
#
#     plt.rc('text', usetex=True)
#     plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
#     plt.rcParams['font.family'] = 'serif'
#     plt.rcParams['font.serif'] = 'Latin Modern Roman'
#     plt.rcParams['xtick.labelsize'] = ticks_font_size
#     plt.rcParams['ytick.labelsize'] = ticks_font_size
#
#     dic_keys = ['sun', 'nun', 'un']
#     bar_labels = ['Sun reference', 'NUN-OR', 'UN-OR']
#
#     fig = plt.figure(figsize=fig_size, dpi=dpi)
#     fig.suptitle(
#         f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$' + f' {round(location.lat, 2)}' + r'$^{\circ})$',
#         fontsize=axes_label_size)
#
#     for i, es in enumerate(sources):
#
#         avg_acceptance_data = avg_acc4designs(geom=geometry, site=location, source=es)
#         ns, ew = avg_acceptance_data
#
#         ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
#         ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)
#
#         for k, dic in enumerate([ns, ew]):
#             ax = ax1 if k == 0 else ax2
#
#             heights = [dic[k] for k in dic_keys]
#
#             heights = array(heights).round(2)
#             bar_pos = arange(len(bar_labels))
#
#             ax.bar(bar_pos, heights, align='center')
#             add_labels(bar_pos, heights,
#                        ax=ax,
#                        font_size=ticks_font_size)
#
#             # Create names on the x-axis
#             ax.set_ylim([0, 1.12 * heights.max()])
#             ax.yaxis.set_ticklabels([])
#
#             if k == 0:
#                 ax.set_ylabel(r'$\bar{\theta}_a$', fontsize=axes_label_size)
#
#             if i == 1:
#                 ax.set_xticks(bar_pos, bar_labels, rotation=0, ha='center', fontsize=axes_label_size)
#             else:
#                 ax.set_xticks([])
#
#             if i == 0:
#                 plot_index = '(a)' if k == 0 else '(b)'
#
#             else:
#                 plot_index = '(c)' if k == 0 else '(d)'
#
#             ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
#                          fontsize=axes_label_size)
#
#     plt.tight_layout()
#     if figure_path is not None:
#         plt.savefig(Path(figure_path, f'{geometry.name}_{location.name}_uniform_avg_acceptances.svg'))
#     plt.show()
#
#     return None
########################################################################################################################
########################################################################################################################


def read_acceptance(geom: lfr_geometry, site: SiteData, source: RadialSource):

    ns_acc, ew_acc = acceptance4designs(geom=geom, site=site, source=source)

    return ns_acc, ew_acc


def add_labels(x, y,
               ax,
               font_size: float):

    for i in range(len(x)):
        ax.text(i, 1.025*y[i], f'{y[i]}' + r"$^{\circ}$", ha='center', fontsize=font_size)

    return None


def plot_non_uniform_acceptances(
        geometry: lfr_geometry, location: SiteData, sources: tuple,
        line_width=1.25,
        fig_width: float = 12, fig_height: float = 12, dpi=450,
        font_size=10,
        figure_path: Path = None, fig_format='pdf'):

    assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'
    lw = abs(line_width)

    seaborn.set(style='whitegrid')
    seaborn.set_context('notebook')

    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Latin Modern Roman'

    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width/2.54, fig_height/2.54))

    fig.suptitle(
        f'LFC {geometry.name[-1]}: {location.name} ('
        + r'$\phi=$' + f' {round(location.lat, 2)}' + r'$^{\circ})$' + '\n',
        fontsize=font_size)

    for i, es in enumerate(sources):

        data = read_acceptance(site=location, geom=geometry, source=es)
        ns, ew = data

        ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
        ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)

        for k, dic in enumerate([ns, ew]):

            ax = ax1 if k == 0 else ax2

            ax.tick_params(labelsize=font_size - 1, pad=-3)

            ax.plot(*dic['sun'].T, color='green', label='Sun reference', lw=lw)
            ax.plot(*dic['spr'].T, color='red', label='Specific reference', lw=lw)
            ax.plot(*dic['boito'].T, color='blue', label='BG design', lw=lw)
            ax.plot(*dic['nun'].T, color='magenta', label='NUN-OR', lw=lw)

            if i == 0:
                plot_index = '(a)' if k == 0 else '(b)'
            else:
                plot_index = '(c)' if k == 0 else '(d)'

            ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
                         fontsize=font_size)

            if k == 0:
                ax.set_ylabel(r"$\theta_{a}$", fontsize=font_size)

            if i == 0:
                ax.xaxis.set_ticklabels([])
            else:
                ax.set_xlabel(r"$\theta_{T}$", fontsize=font_size)
                ax.xaxis.set_major_formatter(EngFormatter(unit=r"$^{\circ}$", sep=""))

            ax.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.2f}$^\circ$"))

        if i == 0:
            ax1.legend(fontsize=font_size - 1)

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'{geometry.name}_{location.name}_non_uniform_acceptances.{fig_format}'))
    plt.show()

    return None


def plot_non_uniform_avg_acceptances(
        geometry: lfr_geometry, location: SiteData, sources: tuple,
        fig_width: float = 12, fig_height: float = 12, dpi=450,
        font_size=8,
        figure_path: Path = None, fig_format='pdf'):

    seaborn.set(style='darkgrid')
    seaborn.set_context('notebook')

    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Latin Modern Roman'

    dic_keys = ['sun', 'spr', 'boito', 'nun']
    bar_labels = ['Sun\nreference', 'Specific\nreference', 'BG design', 'NUN-OR']

    fig = plt.figure(figsize=(fig_width/2.54, fig_height/2.54),
                     dpi=dpi)
    fig.suptitle(
        f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$' + f' {round(location.lat, 2)}'
        + r'$^{\circ})$' + '\n',
        fontsize=font_size)

    for i, es in enumerate(sources):

        avg_acceptance_data = avg_acc4designs(geom=geometry, site=location, source=es)
        ns, ew = avg_acceptance_data

        ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
        ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)

        for k, dic in enumerate([ns, ew]):

            ax = ax1 if k == 0 else ax2

            heights = [dic[k] for k in dic_keys]
            heights = array(heights).round(2)
            bar_pos = arange(len(bar_labels))

            ax.bar(bar_pos, heights, align='center')
            add_labels(bar_pos, heights,
                       ax=ax,
                       font_size=font_size - 1)

            # Create names on the x-axis
            ax.set_ylim([0, 1.12*heights.max()])
            ax.yaxis.set_ticklabels([])

            if k == 0:
                ax.set_ylabel(r'$\bar{\theta}_a$', fontsize=font_size)

            if i == 1:
                ax.set_xticks(bar_pos, bar_labels, rotation=0, ha='center', fontsize=font_size-1)
            else:
                ax.set_xticks([])

            if i == 0:
                plot_index = '(a)' if k == 0 else '(b)'

            else:
                plot_index = '(c)' if k == 0 else '(d)'

            ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
                         fontsize=font_size)

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'{geometry.name}_{location.name}_non_uniform_avg_acceptances.{fig_format}'))
    plt.show()

    return None


def plot_uniform_avg_acceptances(
        geometry: lfr_geometry, location: SiteData, sources: tuple,
        fig_width: float = 10, fig_height: float = 6, dpi=300,
        font_size=8,
        figure_path: Path = None, fig_format='pdf'):

    seaborn.set(style='darkgrid')
    seaborn.set_context('notebook')

    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Latin Modern Roman'

    dic_keys = ['sun', 'nun', 'un']
    bar_labels = ['Sun\nreference', 'NUN-OR', 'UN-OR']

    fig = plt.figure(figsize=(fig_width/2.54, fig_height/2.54), dpi=dpi)
    fig.suptitle(
        f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$' + f' {round(location.lat, 2)}'
        + r'$^{\circ})$' + '\n',
        fontsize=font_size)

    for i, es in enumerate(sources):

        avg_acceptance_data = avg_acc4designs(geom=geometry, site=location, source=es)
        ns, ew = avg_acceptance_data

        ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
        ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)

        for k, dic in enumerate([ns, ew]):
            ax = ax1 if k == 0 else ax2

            ax.tick_params(labelsize=font_size - 1)

            heights = [dic[k] for k in dic_keys]

            heights = array(heights).round(2)
            bar_pos = arange(len(bar_labels))

            ax.bar(bar_pos, heights, align='center')
            add_labels(bar_pos, heights,
                       ax=ax,
                       font_size=font_size)

            # Create names on the x-axis
            ax.set_ylim([0, 1.12 * heights.max()])
            ax.yaxis.set_ticklabels([])

            if k == 0:
                ax.set_ylabel(r'$\bar{\theta}_a$', fontsize=font_size)

            if i == 1:
                ax.set_xticks(bar_pos, bar_labels, rotation=0, ha='center', fontsize=font_size - 1)
            else:
                ax.set_xticks([])

            if i == 0:
                plot_index = '(a)' if k == 0 else '(b)'

            else:
                plot_index = '(c)' if k == 0 else '(d)'

            ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
                         fontsize=font_size)

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'{geometry.name}_{location.name}_uniform_avg_acceptances.{fig_format}'))
    plt.show()

    return None


########################################################################################################################

#################################################################################
#################################################################################


# A folder to save the figures
figures_path = Path(Path.cwd(), 'Figures', 'acceptance')
figures_path.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':

    # plot_non_uniform_acceptances(
    #     geometry=lfc_1, location=evora,
    #     sources=(ES1, ES2),
    #     fig_width=12, fig_height=10,
    #     font_size=8, line_width=0.75,
    #     figure_path=figures_path)

    # plot_non_uniform_avg_acceptances(
    #     geometry=lfc_1, location=evora,
    #     sources=(ES1, ES2),
    #     fig_width=13, fig_height=8,
    #     font_size=8,
    #     figure_path=figures_path)

    # plot_uniform_avg_acceptances(
    #     geometry=lfc_1, location=evora,
    #     fig_width=11, fig_height=8, font_size=8,
    #     sources=(ES1, ES2),
    #     figure_path=figures_path, fig_format='pdf')

    pass
