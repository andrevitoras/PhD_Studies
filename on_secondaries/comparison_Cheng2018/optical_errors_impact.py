import time
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn

from numpy import sqrt, linspace, zeros, absolute, array
from scipy.interpolate import interp1d
from tqdm import tqdm

from solaren.scopy import OpticalProperty
from solaren.scopy.sunlight import RadialSource, SiteData
from solaren.soltracepy import Trace
from utils import dic2json

from on_secondaries.comparison_Cheng2018.Cheng2018_classes_and_functions import lfc_optic, get_optical_analysis_data

from on_secondaries.comparison_Cheng2018.Cheng2018_base_models import (sun_source, slope_error, tracking_error,
                                                                       OpticalSettings, primaries_reflectivity,
                                                                       secondary_reflectivity, tube_absorptivity,
                                                                       glass_cover_transmissivity, cheng_cpc, cheng_cec,
                                                                       cheng_lfc, cheng_loc)

from on_secondaries.comparison_Cheng2018.Cheng2018_simulations import simulation_path, trace_options, dni

########################################################################################################################

seaborn.set(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

ticks_font_size = 14
plt.rcParams['xtick.labelsize'] = ticks_font_size
plt.rcParams['ytick.labelsize'] = ticks_font_size

axes_label_size = 14
plt.rcParams['axes.labelsize'] = axes_label_size
plt.rcParams['axes.titlesize'] = axes_label_size


########################################################################################################################
# Functions for this program ###########################################################################################


def generate_settings(overall_error: float):
    primaries_reflectance = OpticalProperty.reflector(name='primaries_refl',
                                                      rho=primaries_reflectivity,
                                                      slope_error=0., spec_error=overall_error)

    secondary_reflectance = OpticalProperty.secondary(name='secondary_refl',
                                                      rho=secondary_reflectivity,
                                                      slope_error=slope_error, spec_error=0.)

    absorber_properties = OpticalProperty.evacuated_tube(alpha=tube_absorptivity,
                                                         tau=glass_cover_transmissivity,
                                                         ref_index=1.52)

    optical_properties = OpticalSettings(primaries_property=primaries_reflectance,
                                         secondary_property=secondary_reflectance,
                                         absorber_property=absorber_properties)

    return optical_properties


def error_analysis(lfc: lfc_optic, dni_value: float,
                   file_name: str, file_path: Path,
                   source: RadialSource, trace: Trace):

    sigma0 = sqrt(4 * slope_error ** 2 + 4 * tracking_error ** 2)
    errors_values = linspace(start=sigma0, stop=5 * sigma0, num=10)

    results = {}

    now = datetime.now()
    print(f'It is now {now}, and error analysis on flux simulations for {lfc.name} has began.')
    for i, error in enumerate(tqdm(errors_values)):

        name = f'{file_name}_{round(1000 * error, 2)}'
        optics = generate_settings(overall_error=error)

        eta, delta = get_optical_analysis_data(lfc=lfc, dni=dni_value,
                                               file_name=name,
                                               file_path=file_path,
                                               optics=optics,
                                               source=source, trace=trace)

        results[error] = eta

    time.sleep(1.)
    return results


def avg_values(lfc: lfc_optic, dni_value: float,
               file_path: Path,
               source: RadialSource, trace: Trace, site_data: SiteData):

    # the transversal and longitudinal incidence angles for the location
    transversal_angles, longitudinal_angles = site_data.linear_angles(NS=True, solar_longitudinal=False)

    # adding the linear angles to the location tmy
    # concentrator is symmetric
    tmy = site_data.tmy_data
    tmy['theta_t'] = absolute(transversal_angles)
    tmy['theta_l'] = absolute(longitudinal_angles)

    # ns_calculations
    ns_tmy = tmy[tmy['theta_t'] <= 85.]
    ns_angles = ns_tmy['theta_t']
    ns_dni = ns_tmy['dni'].array

    # ew_calculations
    ew_tmy = tmy[tmy['theta_l'] <= 85.]
    ew_angles = ew_tmy['theta_l']
    ew_dni = ew_tmy['dni'].array

    error_data = error_analysis(lfc=lfc, dni_value=dni_value,
                                file_name='optic', file_path=file_path,
                                source=source, trace=trace)

    error_values = array(list(error_data.keys()))
    ns_avg = zeros(shape=(error_values.shape[0], 2))
    ew_avg = zeros(shape=(error_values.shape[0], 2))

    for i, err in enumerate(error_values):
        efficiency_data = error_data[err]
        eff_func = interp1d(*efficiency_data.T, kind='linear')

        ns_avg[i] = err*1000, (eff_func(ns_angles).dot(ns_dni) / ns_dni.sum()).round(5)
        ew_avg[i] = err*1000, (eff_func(ew_angles).dot(ew_dni) / ew_dni.sum()).round(5)

    json_file_path = Path(Path.cwd(), 'error_impact_files')
    json_file_path.mkdir(parents=True, exist_ok=True)

    export_dic = {'ns': ns_avg, 'ew': ew_avg}

    dic2json(d=export_dic, file_path=json_file_path, file_name=f'{lfc.name}')

    return ns_avg, ew_avg


########################################################################################################################
########################################################################################################################

files_path = Path(simulation_path, 'optical_errors_impact')

if __name__ == '__main__':

    optics_list = [cheng_lfc, cheng_cpc, cheng_cec]
    labels = ["Cheng's optic", 'CPC', 'CEC']

    # for lfr in optics_list:
    #     data = error_analysis(lfc=lfr, dni_value=dni,
    #                           file_name='optic', file_path=files_path,
    #                           source=sun_source, trace=trace_options)

    fig = plt.figure(dpi=300, figsize=(10, 4))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    for lfr, label in zip(optics_list, labels):

        ns, ew = avg_values(lfc=lfr, dni_value=dni,
                            file_path=files_path,
                            source=sun_source, trace=trace_options, site_data=cheng_loc)

        ax1.plot(*ns.T, label=label)
        ax2.plot(*ew.T)

    ax1.axes.set_xlabel(r'$\sigma_{o}$ [mrad]')
    ax2.axes.set_xlabel(r'$\sigma_{o}$ [mrad]')

    ax1.axes.set_ylabel(r'$\bar{\eta}$')

    ax1.set_title('(a) NS-mounting')
    ax2.set_title('(b) EW-mounting')

    ax1.legend()
    plt.tight_layout()
    plt.show()




