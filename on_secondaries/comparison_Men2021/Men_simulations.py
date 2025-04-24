
from on_secondaries.comparison_Men2021.Men_base_models import *

simulation_path = Path(Path.cwd(), 'simulation_files')
simulation_path.mkdir(parents=True, exist_ok=True)

optic_number_index = 15

men_optic, \
    men_cpc, men_cpc_red_gap, \
    men_cec, men_cec_red_gap = get_men_optics(index_number=optic_number_index,
                                              secondary_contour_points=contour_points_secondary_optics)

# calculating the number of rays intersections ##############
aperture_area = (dst(men_optic.f1, men_optic.f2) * men_optic.concentrator_length) * 1e-6  # in m2.
rays_density = 9413
ray_intersections = int(rays_density * aperture_area)
ray_intersections = ray_intersections if ray_intersections > 2e5 else 2e5

# print(ray_intersections)
##############################################################

# settings for the trace of rays #############################################
dni = 1000.
trace_options = Trace(rays=ray_intersections, cpus=10, seed=123,
                      sunshape=True, optical_errors=True, simulate=True)
##################################################################################

# for the acceptance simulations
angle_step = 0.05
norm_flux_lower_value = 0.85

########################################################################################################################
# SIMULATIONS ##########################################################################################################

if __name__ == '__main__':

    # Flux simulations ############################################################################

    # Men2021 optic
    _, _ = optical_analysis(lfc=men_optic,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            trace=trace_options, soltrace_version=2012)

    # Men2021 with an edge-rays CPC
    _, _ = optical_analysis(lfc=men_cpc,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            trace=trace_options, soltrace_version=2012)

    _, _ = optical_analysis(lfc=men_cpc_red_gap,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            trace=trace_options, soltrace_version=2012)

    # Men2021 with an edge-rays CEC
    # _, _ = optical_analysis(lfc=men_cec,
    #                         file_name='optic', file_path=simulation_path, dni=dni,
    #                         source=sun_source, optics=optical_properties,
    #                         trace=trace_options, soltrace_version=2012)
    #
    # _, _ = optical_analysis(lfc=men_cec_red_gap,
    #                         file_name='optic', file_path=simulation_path, dni=dni,
    #                         source=sun_source, optics=optical_properties,
    #                         trace=trace_options, soltrace_version=2012)

    ###################################################################################################

    # Acceptance simulations ##########################################################################

    _ = acceptance_analysis(lfc=men_optic,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            dt=angle_step, lower_value=norm_flux_lower_value,
                            trace=trace_options, soltrace_version=2012)

    # Men2021 with an edge-rays CPC
    _ = acceptance_analysis(lfc=men_cpc,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            dt=angle_step, lower_value=norm_flux_lower_value,
                            trace=trace_options, soltrace_version=2012)

    _ = acceptance_analysis(lfc=men_cpc_red_gap,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            dt=angle_step, lower_value=norm_flux_lower_value,
                            trace=trace_options, soltrace_version=2012)

    # Men2021 with an edge-rays CEC
    # _ = acceptance_analysis(lfc=men_cec,
    #                         file_name='optic', file_path=simulation_path, dni=dni,
    #                         source=sun_source, optics=optical_properties,
    #                         dt=angle_step, lower_value=norm_flux_lower_value,
    #                         trace=trace_options, soltrace_version=2012)
    #
    # _ = acceptance_analysis(lfc=men_cec_red_gap,
    #                         file_name='optic', file_path=simulation_path, dni=dni,
    #                         source=sun_source, optics=optical_properties,
    #                         dt=angle_step, lower_value=norm_flux_lower_value,
    #                         trace=trace_options, soltrace_version=2012)

    ###################################################################################################

########################################################################################################################
########################################################################################################################
