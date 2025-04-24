import json
from pathlib import Path

from matplotlib.ticker import EngFormatter, PercentFormatter, StrMethodFormatter
from niopy.geometric_transforms import dst

import seaborn
from matplotlib import pyplot as plt
from numpy import array

from solaren.scopy.sunlight import RadialSource, SiteData
from on_the_primaries_radius.base_models import evora, aswan, lfc_1, lfc_2, ES1, ES2, lfr_geometry

######################################################
######################################################


######################################################
# Plot settings ######################################
seaborn.set_theme(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

########################################################################################################################
# Functions for this program  ##########################################################################################


def read_etas(geom: lfr_geometry, site: SiteData, source: RadialSource):
    file_path = Path(Path.cwd(), site.name, geom.name, source.name, 'etas_data.json')

    with open(file_path, 'r') as opened_file:
        etas_data = json.load(opened_file)

    for k in etas_data['ns'].keys():
        etas_data['ns'][k] = array(etas_data['ns'][k])
        etas_data['ew'][k] = array(etas_data['ew'][k])

    return etas_data

# OLD PLOT FUNCTIONS ###################################################################################################

# def plot_non_uniform_etas(geometry: lfr_geometry, location: SiteData, sources: tuple,
#                           figure_path: Path = None,
#                           line_width=1.5, marker_size=8.,
#                           fig_size=(12, 8), dpi=450,
#                           legend_font_size=14, ticks_font_size=14, axes_label_size=18):
#
#     assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'
#
#     lw = abs(line_width)
#     ms = abs(marker_size)
#
#     plt.rcParams['xtick.labelsize'] = ticks_font_size
#     plt.rcParams['ytick.labelsize'] = ticks_font_size
#
#     fig = plt.figure(figsize=fig_size, dpi=dpi)
#     fig.suptitle(
#         f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$' + f' {round(location.lat, 2)}' + r'$^{\circ})$',
#         fontsize=axes_label_size)
#
#     for i, es in enumerate(sources):
#
#         data = read_etas(site=location, geom=geometry, source=es)
#         ns, ew = data['ns'], data['ew']
#
#         ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
#         ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)
#
#         for k, dic in enumerate([ns, ew]):
#             ax = ax1 if k == 0 else ax2
#             ax.plot(*dic['dp'].T, color='black', label=r"Rabl's design: $\bar{\eta}(\theta_{d})$", lw=lw)
#             ax.plot(*dic['opt_dp'].T, color='orange', label=r'$\bar{\eta}^{*}$', lw=0, marker='s', ms=ms)
#             ax.axhline(y=dic['spr'], color='red', label='Specific reference', lw=lw)
#             ax.plot(*dic['sun'].T, color='green', label='Sun reference', lw=0, marker='.', ms=ms + 4)
#             ax.axhline(y=dic['boito'], color='blue', label='BG design', lw=lw)
#             ax.axhline(y=dic['nun'], color='magenta', label='NUN-OR', lw=lw)
#
#             if i == 0:
#                 plot_index = '(a)' if k == 0 else '(b)'
#             else:
#                 plot_index = '(c)' if k == 0 else '(d)'
#
#             ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
#                          fontsize=legend_font_size + 2)
#             # ax.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.4f}"))
#
#             if k == 0:
#                 ax.set_ylabel(r"$\bar{\eta}$", fontsize=axes_label_size)
#
#             if i == 0:
#                 ax.xaxis.set_ticklabels([])
#             else:
#                 ax.set_xlabel(r"$\theta_{d}$", fontsize=axes_label_size)
#                 ax.xaxis.set_major_formatter(EngFormatter(unit=r"$^{\circ}$", sep=""))
#
#         if i == 0:
#             ax2.legend(fontsize=legend_font_size)
#
#     plt.tight_layout()
#     if figure_path is not None:
#         plt.savefig(Path(figure_path, f'{location.name}_{geometry.name}_non_uniform_etas.svg'))
#     plt.show()
#
#     return None


# def plot_var_non_uniform_etas(geometry: lfr_geometry, location: SiteData, sources: tuple,
#                               fig_size=(12, 5), dpi=300, figure_path: Path = None,
#                               bar_width: float = 0.2,
#                               legend_font_size=14, ticks_font_size=14, axes_label_size=16):
#     assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'
#
#     plt.rcParams['xtick.labelsize'] = ticks_font_size
#     plt.rcParams['ytick.labelsize'] = ticks_font_size
#
#     fig = plt.figure(figsize=fig_size, dpi=dpi)
#     fig.suptitle(f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$'
#                  + f' {round(location.lat, 2)}' + r'$^{\circ})$')
#
#     for i, es in enumerate(sources):
#
#         ax = fig.add_subplot(1, 2, 1 if i == 0 else 2)
#
#         data = read_etas(site=location, geom=geometry, source=es)
#         ns, ew = data['ns'], data['ew']
#
#         ns_data = array([ns['nun'], ns['opt_dp'][0][1], ns['sun'][0][1], ns['spr'], ns['boito']])
#         ew_data = array([ew['nun'], ew['opt_dp'][0][1], ew['sun'][0][1], ew['spr'], ew['boito']])
#
#         ns_ref = ns_data[0]
#         ew_ref = ew_data[0]
#
#         ns_data = 100 * (ns_data - ns_ref) / ns_ref
#         ew_data = 100 * (ew_data - ew_ref) / ew_ref
#
#         x = array([0, 1.5])
#         y1 = [ns_data[1], ew_data[1]]
#         y2 = [ns_data[2], ew_data[2]]
#         y3 = [ns_data[3], ew_data[3]]
#         y4 = [ns_data[4], ew_data[4]]
#
#         ax.bar(x - bar_width, y1, bar_width, label=r'$\theta_{d}^{*}$', color='orange')
#         ax.bar(x, y2, bar_width, label='Sun reference', color='green')
#         ax.bar(x + bar_width, y3, bar_width, label='Specific reference', color='red')
#         ax.bar(x + 2 * bar_width, y4, bar_width, label='BG design', color='blue')
#         ax.set_xticks(x + 0.5 * bar_width, ['NS-mounting', 'EW-mounting'])
#         ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
#         ax.yaxis.set_major_formatter(PercentFormatter(decimals=1))
#         ax.set_title(f'(a) {es.name}' if i == 0 else f'(b) {es.name}', y=-0.1, fontsize=legend_font_size + 2)
#
#         if i == 1:
#             ax.legend(ncol=2, fontsize=legend_font_size)
#         else:
#             ax.set_ylabel(r"$\Delta\bar{\eta}$", fontsize=axes_label_size)
#
#     plt.tight_layout()
#     if figure_path is not None:
#         plt.savefig(Path(figure_path, f'{location.name}_{geometry.name}_non_uniform_etas_variation.svg'))
#     plt.show()
#
#     return None

# def plot_var_uniform_etas(geometry: lfr_geometry, locations: tuple, sources: tuple,
#                           fig_size=(12, 8), dpi=300, figure_path: Path = None,
#                           bar_width: float = 0.1,
#                           legend_font_size=14, ticks_font_size=14, axes_label_size=16):
#
#     assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'
#     assert len(locations) <= 2, 'More than two locations are being inputted. Please, check function arguments!'
#
#     plt.rcParams['xtick.labelsize'] = ticks_font_size
#     plt.rcParams['ytick.labelsize'] = ticks_font_size
#
#     fig = plt.figure(figsize=fig_size, dpi=dpi)
#     fig.suptitle(f'LFC {geometry.name[-1]}', fontsize=axes_label_size+2)
#
#     for i, loc in enumerate(locations):
#
#         ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
#         ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)
#
#         for k, es in enumerate(sources):
#
#             data = read_etas(site=loc, geom=geometry, source=es)
#             ns, ew = data['ns'], data['ew']
#
#             ax = ax1 if k == 0 else ax2
#
#             ns_data = array([ns['nun'], ns['sun'][0][1], ns['un'][1]])
#             ew_data = array([ew['nun'], ew['sun'][0][1], ew['un'][1]])
#
#             ns_ref = ns_data[0]
#             ew_ref = ew_data[0]
#
#             ns_data = 100 * (ns_data - ns_ref) / ns_ref
#             ew_data = 100 * (ew_data - ew_ref) / ew_ref
#
#             x = array([0.5, 1])
#             y1 = [ns_data[1], ew_data[1]]
#             y2 = [ns_data[2], ew_data[2]]
#
#             ax.bar(x - bar_width, y1, bar_width, label='Sun reference', color='green')
#             ax.bar(x, y2, bar_width, label='UN-OR', color='navy')
#
#             ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
#             ax.set_xticks(x - 0.5 * bar_width, ['NS-mounting', 'EW-mounting'], fontsize=ticks_font_size)
#
#             ax.yaxis.set_major_formatter(PercentFormatter(decimals=1))
#
#             if i == 0:
#                 plot_index = '(a)' if k == 0 else '(b)'
#                 if k == 1:
#                     ax.legend(fontsize=legend_font_size)
#
#             else:
#                 plot_index = '(c)' if k == 0 else '(d)'
#
#             ax.set_title(f'{plot_index} {loc.name}: {es.name}', y=-0.1, fontsize=axes_label_size)
#
#             if k == 0:
#                 ax.set_ylabel(r"$\Delta\bar{\eta}$", fontsize=axes_label_size)
#
#             if i == 1:
#                 ax.set_xticks([])
#
#     plt.tight_layout()
#     if figure_path is not None:
#         plt.savefig(Path(figure_path, f'{geometry.name}_uniform_eta_variation.svg'))
#     plt.show()
#
#     return None
#
#
# def plot_uniform_data(geometry: lfr_geometry, location: SiteData, sources: tuple,
#                       fig_size=(12, 7), dpi=300, figure_path: Path = None,
#                       line_width=1.5, marker_size=8.,
#                       legend_font_size=14, ticks_font_size=14, axes_label_size=16):
#
#     assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'
#
#     lw = abs(line_width)
#     ms = abs(marker_size)
#
#     plt.rcParams['xtick.labelsize'] = ticks_font_size
#     plt.rcParams['ytick.labelsize'] = ticks_font_size
#
#     fig = plt.figure(figsize=fig_size, dpi=dpi)
#     fig.suptitle(f'LFC {geometry.name[-1]}: {location.name} ('
#                  + r'$\phi=$' + f' {round(location.lat, 2)}' + r'$^{\circ})$')
#
#     for i, es in enumerate(sources):
#
#         data = read_etas(site=location, geom=geometry, source=es)
#         ns, ew = data['ns'], data['ew']
#
#         ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
#         ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)
#
#         for k, dic in enumerate([ns, ew]):
#             ax = ax1 if k == 0 else ax2
#
#             f_distances = [dst(hc, geometry.rec_aim) for hc in geometry.centers]
#             min_f = min(f_distances)
#             max_f = max(f_distances)
#
#             ax.plot(dic['un_data'].T[0] / 1000, dic['un_data'].T[1], lw=lw)
#             ax.plot(dic['un'].T[0] / 1000, dic['un'].T[1], label='Optimum', lw=0, marker='s', ms=ms, color='orange')
#
#             ax.axvline(x=2 * min_f / 1000, ymax=1, label='Specific reference range', color='red', ls='dashed', lw=lw)
#             ax.axvline(x=2 * max_f / 1000, ymax=1, color='red', ls='dashed', lw=lw)
#
#             if k == 0:
#                 ax.set_ylabel(r"$\bar{\eta}$", fontsize=axes_label_size)
#                 if i == 0:
#                     ax.legend(fontsize=legend_font_size)
#
#             if i == 0:
#                 plot_index = '(a)' if k == 0 else '(b)'
#                 ax.set_xticklabels([])
#             else:
#                 plot_index = '(c)' if k == 0 else '(d)'
#                 ax.set_xlabel('Uniform curvature radius [m]', fontsize=axes_label_size)
#
#             ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
#                          fontsize=legend_font_size + 2)
#
#     plt.tight_layout()
#     if figure_path is not None:
#         plt.savefig(Path(figure_path, f'{location.name}_{geometry.name}_uniform_design_data.svg'))
#     plt.show()
#
#     return None

########################################################################################################################


def plot_non_uniform_etas(
        geometry: lfr_geometry, location: SiteData, sources: tuple,
        fig_width=12., fig_height=12., dpi=300,
        font_size=8., line_width=1., marker_size=4.,
        figure_path: Path = None, fig_format='pdf'):

    assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'
    lw = abs(line_width)
    ms = abs(marker_size)

    fig = plt.figure(figsize=(fig_width/2.54, fig_height/2.54),
                     dpi=dpi)
    fig.suptitle(
        f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$' + f' {round(location.lat, 2)}' + r'$^{\circ})$'
        + '\n',
        fontsize=font_size)

    for i, es in enumerate(sources):
        data = read_etas(site=location, geom=geometry, source=es)
        ns, ew = data['ns'], data['ew']

        ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
        ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)

        for k, dic in enumerate([ns, ew]):
            ax = ax1 if k == 0 else ax2
            ax.plot(*dic['dp'].T, color='black', label=r"Rabl's design: $\bar{\eta}(\theta_{d})$", lw=lw)
            ax.plot(*dic['opt_dp'].T, color='orange', label=r'$\bar{\eta}^{*}$', lw=0, marker='s', ms=ms)
            ax.axhline(y=dic['spr'], color='red', label='Specific reference', lw=lw)
            ax.plot(*dic['sun'].T, color='green', label='Sun reference', lw=0, marker='.', ms=ms + 4)
            ax.axhline(y=dic['boito'], color='blue', label='BG design', lw=lw)
            ax.axhline(y=dic['nun'], color='magenta', label='NUN-OR', lw=lw)

            ax.tick_params(axis='both', which='major', labelsize=font_size-1)

            if i == 0:
                plot_index = '(a)' if k == 0 else '(b)'
            else:
                plot_index = '(c)' if k == 0 else '(d)'

            ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
                         fontsize=font_size)
            ax.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.4f}"))

            if k == 0:
                ax.set_ylabel(r"$\bar{\eta}$", fontsize=font_size)

            if i == 0:
                ax.xaxis.set_ticklabels([])
            else:
                ax.set_xlabel(r"$\theta_{d}$", fontsize=font_size)
                ax.xaxis.set_major_formatter(EngFormatter(unit=r"$^{\circ}$", sep=""))

        if i == 0:
            ax2.legend(fontsize=font_size - 1)

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'{location.name}_{geometry.name}_non_uniform_etas.{fig_format}'))
    plt.show()

    return None


def plot_var_non_uniform_etas(
        geometry: lfr_geometry, location: SiteData, sources: tuple,
        fig_width=12., fig_height=6., dpi=300,
        font_size=8., bar_width=1.,
        legend_columns=1, legend_ax=1, legend_loc=None,
        figure_path: Path = None, fig_format='pdf'):

    assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'

    fig = plt.figure(figsize=(fig_width/2.54, fig_height/2.54),
                     dpi=dpi)

    fig.suptitle(f'LFC {geometry.name[-1]}: {location.name} (' + r'$\phi=$'
                 + f' {round(location.lat, 2)}' + r'$^{\circ})$' + '\n', fontsize=font_size)

    for i, es in enumerate(sources):

        ax = fig.add_subplot(1, 2, i + 1)

        data = read_etas(site=location, geom=geometry, source=es)
        ns, ew = data['ns'], data['ew']

        ns_data = array([ns['nun'], ns['opt_dp'][0][1], ns['sun'][0][1], ns['spr'], ns['boito']])
        ew_data = array([ew['nun'], ew['opt_dp'][0][1], ew['sun'][0][1], ew['spr'], ew['boito']])

        ns_ref = ns_data[0]
        ew_ref = ew_data[0]

        ns_data = 100 * (ns_data - ns_ref) / ns_ref
        ew_data = 100 * (ew_data - ew_ref) / ew_ref

        x = array([0, 6*bar_width])
        y1 = [ns_data[1], ew_data[1]]
        y2 = [ns_data[2], ew_data[2]]
        y3 = [ns_data[3], ew_data[3]]
        y4 = [ns_data[4], ew_data[4]]

        ax.bar(x - bar_width, y1, bar_width, label=r'$\theta_{d}^{*}$', color='orange')
        ax.bar(x, y2, bar_width, label='Sun reference', color='green')
        ax.bar(x + bar_width, y3, bar_width, label='Specific reference', color='red')
        ax.bar(x + 2 * bar_width, y4, bar_width, label='BG design', color='blue')

        ax.set_xticks(x + 0.5 * bar_width, ['NS-mounting', 'EW-mounting'])
        ax.yaxis.set_major_formatter(PercentFormatter(decimals=1))
        ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False,
                       pad=0.,
                       labelsize=font_size-1)

        ax.set_title(f'(a) {es.name}' if i == 0 else f'(b) {es.name}',
                     y=-0.15,
                     fontsize=font_size)

        if i == legend_ax:
            ax.legend(ncol=legend_columns, fontsize=font_size-1, loc='best' if legend_loc is None else legend_loc)

        if i == 0:
            ax.set_ylabel(r"$\Delta\bar{\eta}$", fontsize=font_size)

    plt.tight_layout(pad=0.25, w_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'{location.name}_{geometry.name}_non_uniform_etas_variation.{fig_format}'))
    plt.show()

    return None


def plot_var_uniform_etas(
        geometry: lfr_geometry, locations: tuple, sources: tuple,
        fig_width=10., fig_height=5., dpi=300,
        bar_width: float = 0.1,
        font_size=8,
        figure_path: Path = None, fig_format='pdf'):

    assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'
    assert len(locations) <= 2, 'More than two locations are being inputted. Please, check function arguments!'

    fig = plt.figure(figsize=(fig_width/2.54, fig_height/2.54),
                     dpi=dpi)

    fig.suptitle(f'LFC {geometry.name[-1]}\n', fontsize=font_size)

    for i, loc in enumerate(locations):

        ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
        ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)

        for k, es in enumerate(sources):

            data = read_etas(site=loc, geom=geometry, source=es)
            ns, ew = data['ns'], data['ew']

            ax = ax1 if k == 0 else ax2

            ax.tick_params()

            ns_data = array([ns['nun'], ns['sun'][0][1], ns['un'][1]])
            ew_data = array([ew['nun'], ew['sun'][0][1], ew['un'][1]])

            ns_ref = ns_data[0]
            ew_ref = ew_data[0]

            ns_data = 100 * (ns_data - ns_ref) / ns_ref
            ew_data = 100 * (ew_data - ew_ref) / ew_ref

            x = array([0, 3 * bar_width])
            y1 = [ns_data[1], ew_data[1]]
            y2 = [ns_data[2], ew_data[2]]

            ax.bar(x - bar_width, y1, bar_width, label='Sun reference', color='green')
            ax.bar(x, y2, bar_width, label='UN-OR', color='navy')

            ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False,
                           labelsize=font_size - 1, pad=0)

            ax.set_xticks(x - 0.5 * bar_width, ['NS-mounting', 'EW-mounting'])

            ax.yaxis.set_major_formatter(PercentFormatter(decimals=1))

            if i == 0:
                plot_index = '(a)' if k == 0 else '(b)'
                if k == 1:
                    ax.legend(fontsize=font_size - 1)

            else:
                plot_index = '(c)' if k == 0 else '(d)'

            ax.set_title(f'{plot_index} {loc.name}: {es.name}',
                         y=-0.15,
                         fontsize=font_size)

            if k == 0:
                ax.set_ylabel(r"$\Delta\bar{\eta}$", fontsize=font_size)

            if i == 1:
                ax.set_xticks([])

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'{geometry.name}_uniform_eta_variation.{fig_format}'))
    plt.show()

    return None


def plot_uniform_data(
        geometry: lfr_geometry, location: SiteData, sources: tuple,
        fig_width=12, fig_height=10, dpi=300,
        font_size: float = 8,
        line_width=1.5, marker_size=8.,
        figure_path: Path = None, fig_format='pdf'):

    assert len(sources) <= 2, 'More than two effective sources are being inputted. Please, check function arguments!'

    lw = abs(line_width)
    ms = abs(marker_size)

    fig = plt.figure(figsize=(fig_width/2.54, fig_height/2.54),
                     dpi=dpi)
    fig.suptitle(f'LFC {geometry.name[-1]}: {location.name} ('
                 + r'$\phi=$' + f' {round(location.lat, 2)}' + r'$^{\circ})$' + '\n',
                 fontsize=font_size)

    for i, es in enumerate(sources):

        data = read_etas(site=location, geom=geometry, source=es)
        ns, ew = data['ns'], data['ew']

        ax1 = fig.add_subplot(2, 2, 1 if i == 0 else 3)
        ax2 = fig.add_subplot(2, 2, 2 if i == 0 else 4)

        for k, dic in enumerate([ns, ew]):

            ax = ax1 if k == 0 else ax2

            ax.tick_params(labelsize=font_size - 1)

            f_distances = [dst(hc, geometry.rec_aim) for hc in geometry.centers]
            min_f = min(f_distances)
            max_f = max(f_distances)

            ax.plot(dic['un_data'].T[0] / 1000, dic['un_data'].T[1], lw=lw)
            ax.plot(dic['un'].T[0] / 1000, dic['un'].T[1], label='Optimum', lw=0, marker='s', ms=ms, color='orange')

            ax.axvline(x=2 * min_f / 1000, ymax=1, label='Specific reference', color='red', ls='dashed', lw=lw)
            ax.axvline(x=2 * max_f / 1000, ymax=1, color='red', ls='dashed', lw=lw)

            if k == 0:
                ax.set_ylabel(r"$\bar{\eta}$", fontsize=font_size)
                if i == 0:
                    ax.legend(fontsize=font_size - 1)

            if i == 0:
                plot_index = '(a)' if k == 0 else '(b)'
                ax.set_xticklabels([])
            else:
                plot_index = '(c)' if k == 0 else '(d)'
                ax.set_xlabel('Uniform curvature radius [m]', fontsize=font_size)

            ax.set_title(f'{plot_index} NS-mounting: {es.name}' if k == 0 else f'{plot_index} EW-mounting: {es.name}',
                         fontsize=font_size)

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'{location.name}_{geometry.name}_uniform_design_data.{fig_format}'))
    plt.show()

    return None

########################################################################################################################
########################################################################################################################


# A folder to save the figures
figures_path = Path(Path.cwd(), 'Figures', 'efficiency')
figures_path.mkdir(parents=True, exist_ok=True)

########################################################################################################################
# The Non-uniform designs ##############################################################################################

if __name__ == '__main__':

    pass

    # NON-UNIFORM PLOTS ################################################################################################
    plot_non_uniform_etas(geometry=lfc_1, location=evora,
                          sources=(ES1, ES2),
                          fig_width=12., fig_height=10., dpi=300,
                          font_size=8., line_width=1., marker_size=4.,
                          figure_path=figures_path, fig_format='pdf')

    plot_var_non_uniform_etas(geometry=lfc_1, location=evora,
                              sources=(ES1, ES2),
                              fig_width=12., fig_height=6., dpi=300,
                              font_size=8., bar_width=1.,
                              figure_path=figures_path, fig_format='pdf')

    plot_var_non_uniform_etas(geometry=lfc_1, location=aswan,
                              sources=(ES1, ES2),
                              fig_width=12.5, fig_height=6., dpi=300,
                              font_size=8., bar_width=1.,
                              legend_columns=1, legend_ax=0, legend_loc='lower right',
                              figure_path=figures_path, fig_format='pdf')

    plot_var_non_uniform_etas(geometry=lfc_2, location=aswan,
                              sources=(ES1, ES2),
                              fig_width=12.5, fig_height=6, dpi=300,
                              font_size=8., bar_width=1.,
                              legend_columns=1, legend_ax=0, legend_loc='lower right',
                              figure_path=figures_path, fig_format='pdf')

    ####################################################################################################################
    ####################################################################################################################

    # Uniform configuration plots ######################################################################################
    plot_var_uniform_etas(geometry=lfc_1, locations=(evora, aswan), sources=(ES1, ES2),
                          fig_width=11, fig_height=9,
                          font_size=8,
                          figure_path=figures_path, fig_format='pdf')

    plot_uniform_data(geometry=lfc_1, location=evora,
                      sources=(ES1, ES2),
                      fig_width=13, fig_height=9, dpi=300,
                      font_size=8,
                      line_width=1., marker_size=3.,
                      figure_path=figures_path, fig_format='pdf')

    ####################################################################################################################
    ####################################################################################################################
