from on_secondaries.comparison_Aplanatic.aplanat_base_models import *

simulation_path = Path(Path.cwd(), 'simulation_files')
simulation_path.mkdir(parents=True, exist_ok=True)

# calculating the number of rays intersections ##############
aperture_area = (dst(aplanat_lfc.f1, aplanat_lfc.f2) * aplanat_lfc.concentrator_length) * 1e-6  # in m2.
rays_density = 9413
ray_intersections = rays_density * aperture_area
ray_intersections = ray_intersections if ray_intersections > 2e5 else 2e5
##############################################################

# settings for the trace of rays #############################################
dni = 1000.
trace_options = Trace(rays=int(ray_intersections), cpus=10, seed=123,
                      sunshape=True, optical_errors=False, simulate=True)
##################################################################################

# for the acceptance simulations
angle_step = 0.05
norm_flux_lower_value = 0.85

########################################################################################################################
# SIMULATIONS ##########################################################################################################

if __name__ == '__main__':

    # Flux simulations ############################################################################

    # Aplanat primary field and aplanat secondary optic
    _, _ = optical_analysis(lfc=aplanat_lfc,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            trace=trace_options, soltrace_version=2012)

    # Aplanat primary field with an edge-rays CPC
    _, _ = optical_analysis(lfc=aplanat_lfc_cpc_red_gap,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            trace=trace_options, soltrace_version=2012)

    # Aplanat primary field with an edge-rays CEC
    _, _ = optical_analysis(lfc=aplanat_lfc_cec_red_gap,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            trace=trace_options, soltrace_version=2012)

    ###################################################################################################

    # Acceptance simulations ##########################################################################

    # Aplanat primary field and aplanat secondary optic
    _ = acceptance_analysis(lfc=aplanat_lfc,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            dt=angle_step, lower_value=norm_flux_lower_value,
                            trace=trace_options, soltrace_version=2012)

    # Aplanat primary field with an edge-rays CPC
    _ = acceptance_analysis(lfc=aplanat_lfc_cpc_red_gap,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            dt=angle_step, lower_value=norm_flux_lower_value,
                            trace=trace_options, soltrace_version=2012)

    # Aplanat primary field with an edge-rays CEC
    _ = acceptance_analysis(lfc=aplanat_lfc_cec_red_gap,
                            file_name='optic', file_path=simulation_path, dni=dni,
                            source=sun_source, optics=optical_properties,
                            dt=angle_step, lower_value=norm_flux_lower_value,
                            trace=trace_options, soltrace_version=2012)

    ##################################################################################################

########################################################################################################################
########################################################################################################################
