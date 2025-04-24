from on_optimum_designs.RQ.functions_data_analysis import *


if __name__ == '__main__':

    figures_path = Path(Path.cwd(), 'RQ1_Figures')
    figures_path.mkdir(parents=True, exist_ok=True)

    plot_comparison(number_mirrors=(12, 25),
                    boundaries=threshold_condition,
                    configurations=['non-uniform', 'variable-radius'],
                    gen='genH',
                    figwidth=12, figheight=6,
                    font_size=8,
                    file_path=figures_path, plot_arrows=True, file_format='pdf')

    plot_shape_and_filling_factors(number_mirrors=(12, 25),
                                   boundaries=threshold_condition,
                                   configurations=['non-uniform', 'variable-radius'],
                                   gen='genH',
                                   figwidth=10, figheight=5.5,
                                   font_size=8,
                                   file_path=figures_path, file_format='pdf')

    plot_parameters_correlation(number_of_mirrors=12,
                                boundaries=threshold_condition,
                                configurations=['non-uniform', 'variable-radius'],
                                gen='genH', dni_sum=evora.dni_sum/1000.,
                                figwidth=13, figheight=12,
                                font_size=8,
                                file_path=figures_path, file_format='pdf')
    ####################################################################################################################
