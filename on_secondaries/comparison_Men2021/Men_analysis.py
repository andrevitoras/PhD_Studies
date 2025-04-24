import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter, StrMethodFormatter

from on_secondaries.comparison_Men2021.Men_base_models import *

from on_secondaries.comparison_Men2021.Men_simulations import (simulation_path,
                                                               trace_options, dni, angle_step, norm_flux_lower_value)


# Functions for this program ###########################################################################################

# def plot_metrics_all_cases(index: int,
#                            fig_width: float = 6, fig_height: float = 10, dpi=300,
#                            font_size: float = 8, line_width: float = 1., marker_size=2.,
#                            ticks_pad=-3, title_pad=3.,
#                            figure_path: Path = None, fig_format='pdf'):
#
#     lw = abs(line_width)
#     ms = abs(marker_size)
#
#     optic, cpc, cpc_red_gap, _, _ = get_men_optics(index_number=index,
#                                                    secondary_contour_points=121)
#
#     optics_label = [f"Men's optic \\#{optic.index}", 'CPC']
#
#     fig = plt.figure(dpi=dpi,
#                      figsize=(fig_width/2.54, fig_height/2.54))
#
#     ax_1 = fig.add_subplot(3, 2, 1)
#     ax_2 = fig.add_subplot(3, 2, 3)
#     ax_3 = fig.add_subplot(3, 2, 5)
#
#     for lfr, label in zip([optic, cpc_red_gap], optics_label):
#         efficiency_data, uniformity_data, acceptance_data = get_simulations_data(lfc=lfr,
#                                                                                  file_name='optic',
#                                                                                  file_path=simulation_path,
#                                                                                  dni=dni,
#                                                                                  source=sun_source,
#                                                                                  optics=optical_properties,
#                                                                                  trace=trace_options,
#                                                                                  angle_step=angle_step,
#                                                                                  lower_value=norm_flux_lower_value,
#                                                                                  soltrace_version=2012)
#
#         ax_1.plot(*efficiency_data.T, label=label, marker='s', ms=ms, lw=lw)
#         ax_2.plot(*uniformity_data.T, marker='s', ms=ms, lw=lw)
#         ax_3.plot(*acceptance_data.T, marker='s', ms=ms, lw=lw)
#
#     ax_1.legend(fontsize=font_size-1)
#
#     ax_1.set_title('(a) Designs with different gaps\n')
#
#     ax_1.axes.set_ylabel(r'$\eta$', fontsize=font_size)
#     ax_2.axes.set_ylabel(r'$\delta_q$', fontsize=font_size)
#     ax_3.axes.set_ylabel(r'$\beta$', fontsize=font_size)
#
#     ax_1.set_xticklabels([])
#     ax_2.set_xticklabels([])
#
#     ax_3.axes.set_xlabel(r'$\theta_{T}$', fontsize=font_size)
#
#     ax_3.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
#     ax_3.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.2f}$^\circ$"))
#
#     ax_4 = fig.add_subplot(3, 2, 2)
#     ax_5 = fig.add_subplot(3, 2, 4)
#     ax_6 = fig.add_subplot(3, 2, 6)
#
#     for lfr, label in zip([optic, cpc], optics_label):
#         efficiency_data, uniformity_data, acceptance_data = get_simulations_data(lfc=lfr,
#                                                                                  file_name='optic',
#                                                                                  file_path=simulation_path,
#                                                                                  dni=dni,
#                                                                                  source=sun_source,
#                                                                                  optics=optical_properties,
#                                                                                  trace=trace_options,
#                                                                                  angle_step=angle_step,
#                                                                                  lower_value=norm_flux_lower_value,
#                                                                                  soltrace_version=2012)
#
#         ax_4.plot(*efficiency_data.T, label=label, marker='s', ms=4)
#         ax_5.plot(*uniformity_data.T, marker='s', ms=4)
#         ax_6.plot(*acceptance_data.T, marker='s', ms=4)
#
#     ax_4.legend(fontsize=font_size-1)
#
#     ax_4.set_title('(b) Designs with same gap', fontsize=font_size, pad=title_pad)
#
#     ax_4.axes.set_ylabel(r'$\eta$', fontsize=font_size)
#     ax_5.axes.set_ylabel(r'$\delta_q$', fontsize=font_size)
#     ax_6.axes.set_ylabel(r'$\beta$', fontsize=font_size)
#
#     ax_4.set_xticklabels([])
#     ax_5.set_xticklabels([])
#
#     ax_6.axes.set_xlabel(r'$\theta_{T}$', fontsize=font_size)
#
#     ax_6.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
#     ax_6.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.2f}$^\circ$"))
#
#     [ax.tick_params(labelsize=font_size - 1, pad=ticks_pad) for ax in [ax_1, ax_2, ax_3, ax_4, ax_5, ax_6]]
#
#     plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
#     if figure_path is not None:
#         plt.savefig(Path(figure_path, f'men_{optic.index}_avg_metrics_full.{fig_format}'))
#     plt.show()
#
#     return None


# def plot_metrics(optics_to_plot: list, optics_labels: list, figure):
#     ax_1 = figure.add_subplot(3, 1, 1)
#     ax_2 = figure.add_subplot(3, 1, 2)
#     ax_3 = figure.add_subplot(3, 1, 3)
#
#     for optic, label in zip(optics_to_plot, optics_labels):
#         efficiency_data, uniformity_data, acceptance_data = get_simulations_data(lfc=optic,
#                                                                                  file_name='optic',
#                                                                                  file_path=simulation_path,
#                                                                                  dni=dni,
#                                                                                  source=sun_source,
#                                                                                  optics=optical_properties,
#                                                                                  trace=trace_options,
#                                                                                  angle_step=angle_step,
#                                                                                  lower_value=norm_flux_lower_value,
#                                                                                  soltrace_version=2012)
#
#         ax_1.plot(*efficiency_data.T, label=label, marker='s', ms=4)
#         ax_2.plot(*uniformity_data.T, marker='s', ms=4)
#         ax_3.plot(*acceptance_data.T, marker='s', ms=4)
#
#     ax_1.legend(fontsize=axes_label_size)
#
#     ax_1.set_title('(a) Transversal optical efficiency')
#     ax_2.set_title('(b) Circumferential flux non-uniformity index')
#     ax_3.set_title('(c) Acceptance half-angle')
#
#     ax_1.axes.set_ylabel(r'$\eta$')
#     ax_2.axes.set_ylabel(r'$\delta_q$')
#     ax_3.axes.set_ylabel(r'$\beta$')
#
#     ax_1.set_xticklabels([])
#     ax_2.set_xticklabels([])
#
#     ax_3.axes.set_xlabel(r'$\theta_{T}$')
#
#     # ax_1.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
#     # ax_2.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
#     ax_3.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
#     ax_3.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.2f}$^\circ$"))
#
#     return None

# def plot_averages_var(optics_to_plot: list, optics_labels: list,
#                       site_data: SiteData, figure):
#     ns, ew = get_averages(optics=optics_to_plot, optics_labels=optics_labels, site_data=site_data)
#
#     ns_base = ns[optics_labels[0]].round(3)
#     ew_base = ew[optics_labels[0]].round(3)
#
#     ns_var = {
#         optics_labels[1]: 100 * ((ns[optics_labels[1]] - ns[optics_labels[0]]) / ns[optics_labels[0]]).round(4),
#         optics_labels[2]: 100 * ((ns[optics_labels[2]] - ns[optics_labels[0]]) / ns[optics_labels[0]]).round(4),
#     }
#
#     ew_var = {
#         optics_labels[1]: 100 * ((ew[optics_labels[1]] - ew[optics_labels[0]]) / ew[optics_labels[0]]).round(4),
#         optics_labels[2]: 100 * ((ew[optics_labels[2]] - ew[optics_labels[0]]) / ew[optics_labels[0]]).round(4),
#     }
#
#     ax_1 = figure.add_subplot(1, 2, 1)
#     ax_2 = figure.add_subplot(1, 2, 2)
#
#     metrics = (r"Efficiency ($\bar{\eta}$)", r"Non-uniformity ($\bar{\delta_q}$)", r"Acceptance ($\bar{\beta}$)")
#     colors = ('darkorange', 'green')
#
#     x = arange(len(metrics))  # the label locations
#     width = 0.3  # the width of the bars
#     multiplier = 0
#
#     for i, (attribute, measurement) in enumerate(ns_var.items()):
#         offset = width * multiplier
#         rects = ax_1.bar(x + offset, measurement, width, label=attribute, color=colors[i])
#         ax_1.bar_label(rects, padding=0.5, fmt=StrMethodFormatter(u"{x:.1f}$\%$"))
#         multiplier += 1
#
#     ax_1.set_title('(a) NS-mounting\n'
#                    f'{optics_labels[0]} base values = ({ns_base[0]}, {ns_base[1]}, {ns_base[2]}$^\circ$)')
#     ax_1.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.0f}$\%$"))
#
#     ax_1.set_xticks(x + width / 2, metrics)
#     ax_1.legend()
#
#     multiplier = 0
#     for i, (attribute, measurement) in enumerate(ew_var.items()):
#         offset = width * multiplier
#         rects = ax_2.bar(x + offset, measurement, width, label=attribute, color=colors[i])
#         ax_2.bar_label(rects, padding=0.5, fmt=StrMethodFormatter(u"{x:.1f}$\%$"))
#         multiplier += 1
#
#     ax_2.set_title('(b) EW-mounting\n'
#                    f'{optics_labels[0]} base values = ({ew_base[0]}, {ew_base[1]}, {ew_base[2]}$^\circ$)')
#     ax_2.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.0f}$\%$"))
#
#     ax_2.set_xticks(x + width / 2, metrics)
#
#     return None


def optic_averages(optic, site_data):
    ns_avs, ew_avs = get_averaged_values(lfc=optic, site_data=site_data,
                                         file_name='optic',
                                         file_path=simulation_path,
                                         dni=dni,
                                         source=sun_source,
                                         optics=optical_properties,
                                         trace=trace_options,
                                         angle_step=angle_step,
                                         lower_value=norm_flux_lower_value,
                                         soltrace_version=2012)

    return ns_avs, ew_avs


def get_averages(optics: list, optics_labels: list, site_data):
    ns_averages = {}
    ew_averages = {}

    for optic, label in zip(optics, optics_labels):
        ns, ew = get_averaged_values(lfc=optic, site_data=site_data,
                                     file_name='optic',
                                     file_path=simulation_path,
                                     dni=dni,
                                     source=sun_source,
                                     optics=optical_properties,
                                     trace=trace_options,
                                     angle_step=angle_step,
                                     lower_value=norm_flux_lower_value,
                                     soltrace_version=2012)

        ns_averages[label] = array([ns['avg_eff'], ns['avg_uni'], ns['avg_acc']]).round(4)
        ew_averages[label] = array([ew['avg_eff'], ew['avg_uni'], ew['avg_acc']]).round(4)

    return ns_averages, ew_averages


def plot_metrics(
        index_number: int, same_gap=False,
        fig_width: float = 14, fig_height: float = 6, dpi=300,
        font_size: float = 8, line_width: float = 1., marker_size=2.,
        ticks_pad=-3, title_pad=3.,
        figure_path: Path = None, fig_format='pdf'):

    base_optic, \
        cpc_optic, cpc_red_gap, \
        cec_optic, cec_red_gap = get_men_optics(index_number=index_number,
                                                secondary_contour_points=contour_points_secondary_optics)
    if same_gap:
        optics_to_plot = [base_optic, cpc_optic, cec_optic]
    else:
        optics_to_plot = [base_optic, cpc_red_gap, cec_red_gap]

    optics_labels = [f"Men's optic \\#{base_optic.index}", 'CPC', 'CEC']

    lw = abs(line_width)
    ms = abs(marker_size)

    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width / 2.54, fig_height / 2.54))

    ax_1 = fig.add_subplot(3, 1, 1)
    ax_2 = fig.add_subplot(3, 1, 2)
    ax_3 = fig.add_subplot(3, 1, 3)

    for optic, label in zip(optics_to_plot, optics_labels):
        efficiency_data, uniformity_data, acceptance_data = get_simulations_data(lfc=optic,
                                                                                 file_name='optic',
                                                                                 file_path=simulation_path,
                                                                                 dni=dni,
                                                                                 source=sun_source,
                                                                                 optics=optical_properties,
                                                                                 trace=trace_options,
                                                                                 angle_step=angle_step,
                                                                                 lower_value=norm_flux_lower_value,
                                                                                 soltrace_version=2012)

        ax_1.plot(*efficiency_data.T, label=label,
                  lw=lw,
                  marker='s', ms=ms)
        ax_2.plot(*uniformity_data.T,
                  lw=lw,
                  marker='s', ms=ms)
        ax_3.plot(*acceptance_data.T.round(3),
                  lw=lw,
                  marker='s', ms=ms)

    ax_1.legend(fontsize=font_size - 1)

    ax_1.set_title('(a) Transversal optical efficiency', fontsize=font_size, pad=title_pad)
    ax_2.set_title('(b) Circumferential flux non-uniformity', fontsize=font_size, pad=title_pad)
    ax_3.set_title('(c) Acceptance half-angle', fontsize=font_size, pad=title_pad)

    ax_3.axes.set_xlabel(r'$\theta_{T}$', fontsize=font_size)

    ax_1.axes.set_ylabel(r'$\eta$', fontsize=font_size)
    ax_2.axes.set_ylabel(r'$\delta_q$', fontsize=font_size)
    ax_3.axes.set_ylabel(r'$\beta$', fontsize=font_size)

    ax_1.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
    ax_2.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))

    ax_3.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
    ax_3.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.2f}$^\circ$"))

    [ax.tick_params(labelsize=font_size - 1, pad=ticks_pad) for ax in [ax_1, ax_2, ax_3]]

    ax_1.set_xticklabels([])
    ax_2.set_xticklabels([])

    [ax.grid(lw=0.5 * lw) for ax in [ax_1, ax_2, ax_3]]
    [[spine.set_linewidth(0.5 * lw) for spine in ax.spines.values()] for ax in [ax_1, ax_2, ax_3]]

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        if same_gap:
            plt.savefig(Path(figure_path,
                             f'men_{optics_to_plot[0].index}_same_gap_optical_characterization.{fig_format}'))
        else:
            plt.savefig(Path(figure_path,
                             f'men_{optics_to_plot[0].index}_min_gap_optical_characterization.{fig_format}'))
    plt.show()

    return None


def plot_averages_var(
        index_number: int, same_gap: bool,
        site_data: SiteData,
        fig_width: float = 14, fig_height: float = 6, dpi=300,
        font_size: float = 8, line_width=1.,
        ticks_pad=-3, title_pad=3.,
        figure_path: Path = None, fig_format='pdf'):

    base_optic, \
        cpc_optic, cpc_red_gap, \
        cec_optic, cec_red_gap = get_men_optics(index_number=index_number,
                                                secondary_contour_points=contour_points_secondary_optics)
    if same_gap:
        optics_to_plot = [base_optic, cpc_optic, cec_optic]
    else:
        optics_to_plot = [base_optic, cpc_red_gap, cec_red_gap]

    optics_labels = [f"Men's optic \\#{base_optic.index}", 'CPC', 'CEC']

    ns, ew = get_averages(optics=optics_to_plot, optics_labels=optics_labels, site_data=site_data)

    ns_base = ns[optics_labels[0]].round(3)
    ew_base = ew[optics_labels[0]].round(3)

    ns_var = {
        optics_labels[1]: 100 * ((ns[optics_labels[1]] - ns[optics_labels[0]]) / ns[optics_labels[0]]).round(4),
        optics_labels[2]: 100 * ((ns[optics_labels[2]] - ns[optics_labels[0]]) / ns[optics_labels[0]]).round(4),
    }

    ew_var = {
        optics_labels[1]: 100 * ((ew[optics_labels[1]] - ew[optics_labels[0]]) / ew[optics_labels[0]]).round(4),
        optics_labels[2]: 100 * ((ew[optics_labels[2]] - ew[optics_labels[0]]) / ew[optics_labels[0]]).round(4),
    }

    lw = abs(line_width)
    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width / 2.54, fig_height / 2.54))

    ax_1 = fig.add_subplot(1, 2, 1)
    ax_2 = fig.add_subplot(1, 2, 2)

    metrics = ("Efficiency\n" + r"($\bar{\eta}$)",
               "Non-uniformity\n" + r"($\bar{\delta_q}$)",
               "Acceptance\n" + r"($\bar{\beta}$)")

    colors = ('darkorange', 'green')

    x = arange(len(metrics))  # the label locations
    width = 0.4  # the width of the bars
    multiplier = 0

    for i, (attribute, measurement) in enumerate(ns_var.items()):
        offset = width * multiplier
        rects = ax_1.bar(x + offset, measurement, width, label=attribute, color=colors[i])
        ax_1.bar_label(rects, padding=0.5, fontsize=font_size - 2,
                       fmt=StrMethodFormatter(u"{x:.1f}$\%$"))
        multiplier += 1

    ax_1.set_title('(a) NS-mounting\n'
                   f'{optics_labels[0]} values = ({ns_base[0]}, {ns_base[1]}, {ns_base[2]}$^\circ$)',
                   fontsize=font_size, pad=title_pad)
    ax_1.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.0f}$\%$"))

    ax_1.set_xticks(x + width / 2, metrics)
    ax_1.legend(fontsize=font_size - 2)
    ax_1.tick_params(labelsize=font_size - 1, pad=ticks_pad)

    ax_1.grid(lw=0.5 * lw)
    [spine.set_linewidth(0.5 * lw) for spine in ax_1.spines.values()]

    multiplier = 0
    for i, (attribute, measurement) in enumerate(ew_var.items()):
        offset = width * multiplier
        rects = ax_2.bar(x + offset, measurement, width, label=attribute, color=colors[i])
        ax_2.bar_label(rects, padding=0.5, fontsize=font_size - 2,
                       fmt=StrMethodFormatter(u"{x:.1f}$\%$"))
        multiplier += 1

    ax_2.set_title('(b) EW-mounting\n'
                   f'{optics_labels[0]} values = ({ew_base[0]}, {ew_base[1]}, {ew_base[2]}$^\circ$)',
                   fontsize=font_size, pad=title_pad)
    ax_2.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.0f}$\%$"))

    ax_2.set_xticks(x + width / 2, metrics)
    ax_2.tick_params(labelsize=font_size - 1, pad=ticks_pad)

    ax_2.grid(lw=0.5 * lw)
    [spine.set_linewidth(0.5 * lw) for spine in ax_2.spines.values()]

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        if same_gap:
            plt.savefig(Path(figure_path,
                             f'men_{optics_to_plot[0].index}_same_gap_avg_metrics.{fig_format}'))
        else:
            plt.savefig(Path(figure_path,
                             f'men_{optics_to_plot[0].index}_min_gap_avg_metrics.{fig_format}'))
    plt.show()

    return None


def plot_metrics_all_cases(
        index_number: int,
        fig_width: float = 14, fig_height: float = 6, dpi=300,
        font_size: float = 8, line_width=1., marker_size=2,
        ticks_pad=-3, title_pad=3.,
        figure_path: Path = None, fig_format='pdf'):

    lw = abs(line_width)
    ms = abs(marker_size)

    optic, cpc, cpc_red_gap, _, _ = get_men_optics(index_number=index_number,
                                                   secondary_contour_points=121)

    optics_label = [f"Men's optic \\#{optic.index}", 'CPC']

    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width/2.54, fig_height/2.54))

    ax_1 = fig.add_subplot(3, 2, 1)
    ax_2 = fig.add_subplot(3, 2, 3)
    ax_3 = fig.add_subplot(3, 2, 5)

    for lfr, label in zip([optic, cpc_red_gap], optics_label):
        efficiency_data, uniformity_data, acceptance_data = get_simulations_data(lfc=lfr,
                                                                                 file_name='optic',
                                                                                 file_path=simulation_path,
                                                                                 dni=dni,
                                                                                 source=sun_source,
                                                                                 optics=optical_properties,
                                                                                 trace=trace_options,
                                                                                 angle_step=angle_step,
                                                                                 lower_value=norm_flux_lower_value,
                                                                                 soltrace_version=2012)

        ax_1.plot(*efficiency_data.T, label=label, marker='s', ms=ms, lw=lw)
        ax_2.plot(*uniformity_data.T, marker='s', ms=ms, lw=lw)
        ax_3.plot(*acceptance_data.T, marker='s', ms=ms, lw=lw)

    ax_1.legend(fontsize=font_size-1)

    ax_1.set_title('(a) Designs with different gaps', fontsize=font_size, pad=title_pad)

    ax_1.axes.set_ylabel(r'$\eta$', fontsize=font_size)
    ax_2.axes.set_ylabel(r'$\delta_q$', fontsize=font_size)
    ax_3.axes.set_ylabel(r'$\beta$', fontsize=font_size)

    ax_1.set_xticklabels([])
    ax_2.set_xticklabels([])

    ax_3.axes.set_xlabel(r'$\theta_{T}$', fontsize=font_size)

    ax_3.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
    ax_3.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.2f}$^\circ$"))

    ax_4 = fig.add_subplot(3, 2, 2)
    ax_5 = fig.add_subplot(3, 2, 4)
    ax_6 = fig.add_subplot(3, 2, 6)

    for lfr, label in zip([optic, cpc], optics_label):
        efficiency_data, uniformity_data, acceptance_data = get_simulations_data(lfc=lfr,
                                                                                 file_name='optic',
                                                                                 file_path=simulation_path,
                                                                                 dni=dni,
                                                                                 source=sun_source,
                                                                                 optics=optical_properties,
                                                                                 trace=trace_options,
                                                                                 angle_step=angle_step,
                                                                                 lower_value=norm_flux_lower_value,
                                                                                 soltrace_version=2012)

        ax_4.plot(*efficiency_data.T, label=label, marker='s', ms=ms, lw=lw)
        ax_5.plot(*uniformity_data.T, marker='s', ms=ms, lw=lw)
        ax_6.plot(*acceptance_data.T, marker='s', ms=ms, lw=lw)

    ax_4.legend(fontsize=font_size-1)

    ax_4.set_title('(b) Designs with same gap', fontsize=font_size, pad=title_pad)

    # ax_4.axes.set_ylabel(fontsize=font_size)
    # ax_5.axes.set_ylabel(fontsize=font_size)
    # ax_6.axes.set_ylabel(fontsize=font_size)

    ax_4.set_xticklabels([])
    ax_5.set_xticklabels([])

    ax_6.axes.set_xlabel(r'$\theta_{T}$', fontsize=font_size)

    ax_6.xaxis.set_major_formatter(EngFormatter(unit=u"$^\circ$", sep=""))
    ax_6.yaxis.set_major_formatter(StrMethodFormatter(u"{x:.2f}$^\circ$"))

    [ax.tick_params(labelsize=font_size - 1, pad=ticks_pad) for ax in [ax_1, ax_2, ax_3, ax_4, ax_5, ax_6]]

    plt.tight_layout(pad=0.25, w_pad=0.5, h_pad=0.5)
    if figure_path is not None:
        plt.savefig(Path(figure_path, f'men_{optic.index}_optical_characterization_full.{fig_format}'))
    plt.show()

    return None


########################################################################################################################
########################################################################################################################


if __name__ == '__main__':
    pass

    emsp_lat, emsp_long = 38.53, -8.0
    evora = SiteData(name='Evora', latitude=emsp_lat, longitude=emsp_long)

    # plot_metrics(
    #     index_number=15, same_gap=False,
    #     fig_width=6, fig_height=10, font_size=7,
    #     line_width=0.75, marker_size=1,
    #     figure_path=figures_path)

    # plot_averages_var(
    #     index_number=15, same_gap=False,
    #     site_data=evora,
    #     fig_width=14, fig_height=6, font_size=8,
    #     line_width=0.75,
    #     figure_path=figures_path)

    # plot_averages_var(
    #     index_number=15, same_gap=True,
    #     site_data=evora,
    #     fig_width=14, fig_height=6, font_size=8,
    #     line_width=0.75,
    #     figure_path=figures_path)

    plot_metrics_all_cases(
        index_number=4,
        fig_width=10, fig_height=10*(7/9),
        font_size=7, line_width=0.5, marker_size=1,
        title_pad=4, ticks_pad=-4,
        figure_path=figures_path
    )

########################################################################################################################



