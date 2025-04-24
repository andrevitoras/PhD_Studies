from on_optimum_designs.RQ.functions_data_analysis import *

if __name__ == '__main__':
    figures_path = Path(Path.cwd(), 'RQ2_Figures')
    figures_path.mkdir(parents=True, exist_ok=True)

    # non-uniform vs. nun-sun-ref ####################################################################################
    plot_comparison(
        number_mirrors=(12, 25),
        boundaries=threshold_condition,
        configurations=['non-uniform', 'nun-sun-ref'],
        gen='genH',
        figwidth=12, figheight=6,
        font_size=8,
        file_name='NUN_curvatures_configurations_comparison',
        file_path=figures_path, plot_arrows=True, file_format='pdf')

    plot_variables_frequency_distribution(
        number_of_mirrors=12,
        configurations=['non-uniform', 'nun-sun-ref'],
        boundaries=threshold_condition,
        gen='genH',
        figwidth=9, figheight=9, font_size=8,
        num_bins=20,
        plot_bounds_lim=True, legend_ax=2,
        file_name='NUN_curvatures_dec_var_freq_distribution',
        file_path=figures_path, file_format='pdf')

    ####################################################################################################################
    # uniform vs. un-ipa-ref
    plot_comparison(
        number_mirrors=(12, 25),
        boundaries=threshold_condition,
        configurations=['uniform', 'un-ipa-ref'],
        gen='genH',
        figwidth=12, figheight=6,
        font_size=8,
        file_name='UN_curvatures_configurations_comparison',
        file_path=figures_path, plot_arrows=True, file_format='pdf')

    plot_variables_frequency_distribution(
        number_of_mirrors=25,
        configurations=['uniform', 'un-ipa-ref'],
        boundaries=threshold_condition,
        gen='genH',
        figwidth=9, figheight=9, font_size=8,
        num_bins=20,
        plot_bounds_lim=True, legend_ax=2,
        file_name='UN_curvatures_dec_var_freq_distribution',
        file_path=figures_path, file_format='pdf')
    ####################################################################################################################
