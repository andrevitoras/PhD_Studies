from numpy import linspace

from an_optical_method.classes_and_functions import *
from an_optical_method.base_models import lfc_1, ES3

convergence_files_path = Path(Path.cwd(), 'convergence_study_files')


# Functions for this program ###########################################################################################

def convergence_of_points(theta_t: float, theta_l: float,
                          base_lfc: uniform_geometry, ES: EffectiveSource,
                          file_path: Path, force_sim: bool):

    nbr_pts_values = linspace(start=11, stop=701, num=50)

    results = zeros(shape=(nbr_pts_values.shape[0], 2))
    for i, n in enumerate(tqdm(nbr_pts_values)):

        lfc = uniform_geometry(nbr_points=int(n),
                               name=base_lfc.name,
                               number_of_mirrors=base_lfc.number_primaries,
                               mirror_width=base_lfc.mirror_width, mirror_shift=base_lfc.mirror_shift,
                               length=base_lfc.length,
                               receiver_height=base_lfc.receiver_height, radius_design=base_lfc.radius_design)

        df = lfc.analytical_evaluation(ES=ES, symmetric=True, force_sim=force_sim,
                                       path=file_path, file_name=f'{int(n)}_points')

        eta = df[(df['theta_t'] == theta_t) & (df['theta_l'] == theta_l)]['eta'].values[0]

        results[i] = int(n), eta

    return results


def ray_hits_confidence_interval(rayhits: int, number_of_simulations: int,
                                 theta_t: float, theta_l: float,
                                 lfc: uniform_geometry, base_ES: EffectiveSource,
                                 file_path: Path, force_sim: bool, engine='API'):

    files_path = Path(file_path, lfc.name, base_ES.name, f'{engine}')
    files_path.mkdir(parents=True, exist_ok=True)

    if engine == 'GUI':
        n_cpus = 10
    else:
        if rayhits < 1e6:
            n_cpus = 1
        elif 1e6 <= rayhits < 2.5e6:
            n_cpus = 5
        else:
            n_cpus = 10

    sun = SunSettings(sunshape=base_ES.sun.profile,
                      size=base_ES.sun.sunshape_size,
                      rays=int(rayhits),
                      seed=-1, cpus=n_cpus)

    ES = EffectiveSource(name=base_ES.name,
                         sun=sun,
                         optical=base_ES.optical)

    results = zeros(number_of_simulations)
    for i in range(number_of_simulations):

        eta = raytracing_efficiency(file_name=f'optic_rayhits_{int(rayhits)}_run_{i+1}',
                                    files_path=files_path,
                                    lfc=lfc, ES=ES,
                                    theta_t=theta_t, theta_l=theta_l,
                                    force_sim=force_sim, engine=engine, full_optic=False)

        results[i] = eta

    return results


def convergence_of_ray_hits(theta_t: float, theta_l: float,
                            lfc: uniform_geometry, base_ES: EffectiveSource,
                            file_path: Path, force_sim: bool, engine='API'):

    ray_hit_values = [50e3, 100e3, 250e3, 500e3, 1e6, 1.5e6, 3.0e6, 5.0e6]
    simulations_per_rayhit = 100

    conv_analysis_data = {}

    for nrays in tqdm(ray_hit_values):

        etas = ray_hits_confidence_interval(rayhits=int(nrays), number_of_simulations=simulations_per_rayhit,
                                            theta_t=float(theta_t), theta_l=float(theta_l),
                                            lfc=lfc, base_ES=base_ES,
                                            file_path=file_path, force_sim=force_sim, engine=engine)

        conv_analysis_data[nrays] = etas

    return conv_analysis_data


########################################################################################################################


if __name__ == '__main__':

    convergence_of_points(theta_t=0, theta_l=0,
                          base_lfc=lfc_1, ES=ES3,
                          file_path=convergence_files_path, force_sim=False)

    convergence_of_ray_hits(engine='API',
                            theta_t=0., theta_l=0.,
                            lfc=lfc_1, base_ES=ES3,
                            file_path=convergence_files_path, force_sim=False)

    convergence_of_ray_hits(engine='API',
                            theta_t=85., theta_l=85.,
                            lfc=lfc_1, base_ES=ES3,
                            file_path=convergence_files_path, force_sim=False)
