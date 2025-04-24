from numpy import linspace
from utils import dic2json

from on_secondaries.comparison_Men2021.Men_base_models import *
from on_secondaries.comparison_Men2021.Men_simulations import trace_options


def get_optical_analysis_data(file_name: str, file_path: Path,
                              dni: float,
                              lfc: lfc_optic,
                              source: RadialSource,
                              optics: OpticalSettings,
                              trace: Trace,
                              symmetric=False,
                              rec_dx=0., rec_dy=0.,
                              soltrace_version=2012):

    files_path = Path(file_path, 'flux', lfc.name)
    files_path.mkdir(parents=True, exist_ok=True)

    if symmetric:
        angle_range = arange(start=0., stop=90., step=5)
    else:
        angle_range = arange(start=-85, stop=90., step=5)

    mirror_area = (lfc.primary_field.widths.sum() * lfc.concentrator_length) / 1e6  # im m2

    efficiency_data = zeros(shape=(angle_range.shape[0], 2))
    uniformity_data = zeros(shape=(angle_range.shape[0], 2))

    for i, theta in enumerate(angle_range):
        absorber_flux_file_path = simulate_flux(file_name=file_name, file_path=files_path,
                                                theta_t=theta, dni=dni,
                                                lfc=lfc, rec_dx=rec_dx, rec_dy=rec_dy,
                                                source=source, optics=optics, trace=trace,
                                                soltrace_version=soltrace_version)[1]

        if not absorber_flux_file_path.is_file():
            raise ValueError('Absorber flux stats does not exist. Please, run simulation before read stats file!')

        absorber_flux_data = ElementFlux(absorber_flux_file_path)

        efficiency_data[i] = theta, absorber_flux_data.total_flux * 1e3 / (mirror_area * dni)
        uniformity_data[i] = theta, absorber_flux_data.flux_x_uniformity

    return efficiency_data, uniformity_data


def receiver_position_error_analysis(lfc: lfc_optic,
                                     site_data: SiteData,
                                     source: RadialSource,
                                     optics: OpticalSettings,
                                     trace: Trace,
                                     file_path: Path,
                                     symmetric: bool,
                                     rec_dx=0.0, rec_dy=0.):

    efficiency_data, _ = get_optical_analysis_data(file_name=f'optic_{rec_dx}_{rec_dy}',
                                                   file_path=file_path,
                                                   lfc=lfc, symmetric=symmetric,
                                                   rec_dx=rec_dx, rec_dy=rec_dy,
                                                   source=source, optics=optics, trace=trace,
                                                   dni=1000.)

    eff_function = interp1d(*efficiency_data.T, kind='linear')

    # the transversal and longitudinal incidence angles for the location
    transversal_angles, longitudinal_angles = site_data.linear_angles(NS=True, solar_longitudinal=False)

    # adding the linear angles to the location tmy
    # concentrator is symmetric
    tmy = site_data.tmy_data
    # tmy['theta_l'] = longitudinal_angles

    # ns_calculations
    if symmetric:
        tmy['theta_t'] = absolute(transversal_angles)
        ns_tmy = tmy[tmy['theta_t'] <= 85.]
        # angles = ns_tmy['theta_t']
        # dni = ns_tmy['dni'].array
        # ns_eff = eff_function(angles).dot(dni) / dni.sum()
    else:
        tmy['theta_t'] = transversal_angles
        ns_tmy = tmy[absolute(tmy['theta_t']) <= 85.]

    angles = ns_tmy['theta_t']
    dni = ns_tmy['dni'].array

    ns_eff = eff_function(angles).dot(dni) / dni.sum()

    # # ew_calculations
    # ew_tmy = tmy[tmy['theta_l'] <= 85.]
    # angles = ew_tmy['theta_l']
    # dni = ew_tmy['dni'].array
    #
    # ew_eff = eff_function(angles).dot(dni) / dni.sum()

    return ns_eff


def receiver_position_error_range_analysis(lfc: lfc_optic, site_data: SiteData,
                                           source: RadialSource,
                                           optics: OpticalSettings,
                                           trace: Trace,
                                           file_path: Path):

    displacements = linspace(start=-20, stop=20, num=11) * 10

    # displacement_values = array(list(set(displacements.tolist() + [-250., 250.])))
    # displacement_values.sort()

    displacement_values = displacements

    now = datetime.now()
    print(f'It is now {now}, and receiver position error analysis in the x-axis simulations for {lfc.name} has began.')
    avg_eff_by_dx = array([[ds, receiver_position_error_analysis(lfc=lfc, symmetric=False,
                                                                 site_data=site_data,
                                                                 source=source, optics=optics,
                                                                 trace=trace, file_path=file_path,
                                                                 rec_dx=ds, rec_dy=0.)]
                           for ds in tqdm(displacement_values)])
    time.sleep(1.)

    now = datetime.now()
    print(f'It is now {now}, and receiver position error analysis in y-axis simulations for {lfc.name} has began.')
    avg_eff_by_dy = array([[ds, receiver_position_error_analysis(lfc=lfc, symmetric=True,
                                                                 site_data=site_data,
                                                                 source=source, optics=optics,
                                                                 trace=trace, file_path=file_path,
                                                                 rec_dx=0., rec_dy=ds)]
                           for ds in tqdm(displacement_values)])

    ns_dic = {'dx': avg_eff_by_dx, 'dy': avg_eff_by_dy}
    dic2json(d=ns_dic, file_path=file_path, file_name=f'{lfc.name}')

    return avg_eff_by_dx, avg_eff_by_dy

########################################################################################################################
########################################################################################################################


position_error_path = Path(Path.cwd(), 'position_error_files')

emsp_lat, emsp_long = 38.53, -8.0
evora = SiteData(name='Evora', latitude=emsp_lat, longitude=emsp_long)


# lfc_optics = [men_optic, men_cpc_red_gap, men_cec_red_gap]
lfc_optics = [men_optic, men_cpc_red_gap, men_cec_red_gap]

for lfc_concentrator in lfc_optics:
    dx_data, dy_data = receiver_position_error_range_analysis(lfc=lfc_concentrator, site_data=evora,
                                                              file_path=position_error_path,
                                                              source=sun_source, optics=optical_properties,
                                                              trace=trace_options)

    fig = plt.figure(dpi=300)
    ax = fig.add_subplot()

    ax.plot(*dx_data.T, label='X-axis')
    ax.plot(*dy_data.T, label='Y-axis')

    ax.legend()

    ax.set_ylabel(r'$\bar{\eta}$')
    ax.set_xlabel('Receiver displacement [mm]')
    ax.set_title(f'{lfc_concentrator.name}')

    plt.tight_layout()
    plt.show()
