import json
import time
from pathlib import Path

from utils import dic2json

from on_optimum_designs.RQ.classes_and_functions import *

from tqdm import tqdm
from datetime import datetime

from deap import base, creator, tools
from multiprocessing import Pool, cpu_count

########################################################################################################################
# Functions for this program ###########################################################################################

# Fitness and Individual classes
creator.create("FitnessMaxMin", base.Fitness, weights=(1.0, -1.0))
creator.create("Individual", list, fitness=creator.FitnessMaxMin)
############################################################################

# determining the number of cpus in the computer to the parallelization
n_cpus = cpu_count()


class EvolutionaryOperators:

    def __init__(self, n_gens: int, cxpb: float, mupb: float, mustd: float):
        self.n_gens = int(n_gens)
        self.cxpb = abs(cxpb)
        self.mupb = abs(mupb)
        self.mustd = abs(mustd)

        assert 0 < self.cxpb <= 1., 'Crossover probability must be (0,1].'
        assert 0 < self.mupb <= 1., 'Mutation probability must be [0,1).'


def get_lfc_solution(individual: list,
                     number_of_mirrors: int,
                     configuration: str,
                     bounds: VariableBounds):

    decision_vector = [number_of_mirrors] + individual

    if configuration == 'uniform':
        lfc_solution = uniform_geometry(x=decision_vector, bounds=bounds)
    elif configuration == 'non-uniform':
        lfc_solution = non_uniform_geometry(x=decision_vector, bounds=bounds)
    elif configuration == 'variable-radius':
        lfc_solution = variable_radius_geometry(x=decision_vector, bounds=bounds)
    elif configuration == 'nun-sun-ref':
        lfc_solution = nun_sun_ref_geometry(x=decision_vector, bounds=bounds)
    elif configuration == 'un-ipa-ref':
        lfc_solution = un_ipa_ref_geometry(x=decision_vector, bounds=bounds)
    else:
        raise ValueError('Please, select one a possible primary field configuration! See documentation.')

    return lfc_solution


def individual_evaluation(individual: list,
                          number_of_mirrors: int,
                          configuration: str,
                          bounds: VariableBounds,
                          location: SiteData,
                          cum_eff: EffectiveSource,
                          flux_threshold: float):
    # get lfc geometry
    lfc_solution = get_lfc_solution(individual=individual, number_of_mirrors=number_of_mirrors,
                                    configuration=configuration, bounds=bounds)

    # energy collection factor
    ecf = lfc_solution.energy_collection_factor(cum_eff=cum_eff,
                                                location=location,
                                                flux_threshold=flux_threshold)
    # specific cost
    sc = lfc_solution.specific_cost()

    return ecf, sc


def assignRank(individual, rank):
    individual.rank = rank
    return None


def calculate_hypervolume(pareto_individuals: list, reference_point=(-1.0, 0.)) -> float:
    """
        Calculate the hypervolume for a two-objective optimization problem.

        The first objective (efficiency) is to be maximized, and the second objective (cost) is to be minimized.

        Parameters:
        -----------
        solutions : list of individuals objects having the following attribute: individual.fitness.values
            A list of non-dominated Pareto Front solutions, where fitness values are represented as a tuple (f1, f2).
            - f1: Objective 1 value (efficiency, to be maximized)
            - f2: Objective 2 value (cost, to be minimized)
        reference_point : tuple, optional
            The utopian point that dominates all solutions in the objective space.
            Default is (1, 0), representing 100% efficiency and 0 cost.

        Returns:
        --------
        float
            The hypervolume value, which represents the size of the dominated area in the objective space.

        Notes:
        ------
        - Ensure the solutions are non-dominated before calling this function.
        - The solutions should be sorted in descending order of the first objective (f1) and ascending order of the
        second objective (f2).
        """

    # Getting pareto front individuals fitness values
    pareto_front_fitness = [ind.fitness.values for ind in pareto_individuals]

    # Sort Pareto front by decreasing f1 (maximize objective 1)
    sorted_front = sorted(pareto_front_fitness, key=lambda x: -x[0])

    hypervolume = 0.0
    prev_f2 = reference_point[1]  # Start with reference point's f2

    for f1, f2 in sorted_front:
        width = reference_point[0] - f1  # Distance in f1 (maximized)
        height = prev_f2 - f2  # Distance in f2 (minimized)

        # Accumulate the area of the rectangle
        hypervolume += width * height

        # Update the previous f2 value
        prev_f2 = f2

    return hypervolume


def export_population(population: list,
                      file_path: Path,
                      file_name: str):
    pareto_front = [ind for ind in population if ind.rank == 1]
    hv = calculate_hypervolume(pareto_individuals=pareto_front)

    pop_dic = {'population': {'individuals': population,
                              'f1': [ind.fitness.values[0] for ind in population],
                              'f2': [ind.fitness.values[1] for ind in population]},

               'pareto_front': {'individuals': pareto_front,
                                'f1': [ind.fitness.values[0] for ind in pareto_front],
                                'f2': [ind.fitness.values[1] for ind in pareto_front],
                                'hv': hv}
               }

    # exporting population and
    dic2json(d=pop_dic,
             file_path=file_path,
             file_name=file_name)

    return None


def selNSGA2(individuals, k, nd='standard'):
    """
    Modified version of the tools.selNSGA2 function of the DEAP module.
    It is modified to assign a Rank for each pareto front. The Rank = 1 refers to the true Pareto Front.
    """
    if nd == 'standard':
        pareto_fronts = tools.sortNondominated(individuals, k)
    elif nd == 'log':
        pareto_fronts = tools.sortLogNondominated(individuals, k)
    else:
        raise Exception('selNSGA2: The choice of non-dominated sorting '
                        'method "{0}" is invalid.'.format(nd))

    for i, front in enumerate(pareto_fronts):
        [assignRank(individual=ind, rank=i + 1) for ind in front]
        tools.emo.assignCrowdingDist(front)

    chosen = list(tools.emo.chain(*pareto_fronts[:-1]))
    k = k - len(chosen)
    if k > 0:
        sorted_front = sorted(pareto_fronts[-1], key=tools.emo.attrgetter("fitness.crowding_dist"), reverse=True)
        chosen.extend(sorted_front[:k])

    return chosen


def simpleNSGA2(pop,
                n_gens: int,
                toolbox: base.Toolbox,
                cxpb: float, mupb: float, files_path: Path = None) -> list:

    # Defining the starting population: random of user-defined
    if isinstance(pop, (int, float)):
        pop_size = int(pop)
        # Generates a random population of individuals
        population = toolbox.random_population(n=pop_size)
    elif isinstance(pop, (list, ndarray)):
        pop_size = len(pop)
        # pre-defined initial population
        population = toolbox.starting_population(pop)
    else:
        raise ValueError('Please, select a option for the initial population of individuals')

    # Creates the 'select' NSGA-2 evolutionary operator (or overwrite the existing one, if it is already registered)
    if pop_size < 100:
        toolbox.register("select", selNSGA2)
    else:
        toolbox.register("select", selNSGA2, nd='log')

    initial_population_file = Path(files_path, 'gen_0.json')
    if initial_population_file.is_file():
        population = get_individuals_from_file(full_file_path=initial_population_file,
                                               return_pareto=False)
    else:
        # Evaluate initial population
        print(f'Evaluating initial population. {pop_size} individuals')
        with Pool(n_cpus - 2) as pool:
            fitnesses = pool.map(toolbox.evaluate, population)
        for ind, fit in zip(population, fitnesses):
            ind.fitness.values = fit

    # Assign Pareto ranks and crowding distances to the first population
    population = toolbox.select(population, len(population))
    export_population(population=population,
                      file_path=files_path,
                      file_name='gen_0')

    # Loop for the evolutionary process
    print(f'It is {datetime.now()}! Initiating the evolutionary loop.')
    for gen in tqdm(range(n_gens)):
        pop_file_name = f'gen_{gen + 1}'
        population_file = Path(files_path, f'{pop_file_name}.json')

        if population_file.is_file():
            population = get_individuals_from_file(full_file_path=initial_population_file,
                                                   return_pareto=False)
            # Assign Pareto ranks and crowding distances
            population = toolbox.select(population, len(population))
        else:
            # Select parents for mating
            offspring = tools.selTournamentDCD(population, len(population))
            offspring = list(map(toolbox.clone, offspring))

            # Apply crossover and mutation
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < cxpb:
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values
            for mutant in offspring:
                if random.random() < mupb:
                    toolbox.mutate(mutant)
                    del mutant.fitness.values

            # Evaluate offspring with invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            with Pool(n_cpus - 2) as pool:
                fitnesses = pool.map(toolbox.evaluate, invalid_ind)

            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            # Combine parent and offspring populations into a pool of individuals
            individuals_pool = population + offspring

            # Selecting individuals from the pool to the next generation
            population = toolbox.select(individuals_pool, pop_size)
            export_population(population=population,
                              file_path=files_path,
                              file_name=pop_file_name)

    return population


def get_files_path(number_of_mirrors: int,
                   configuration: str,
                   boundary_conditions: BoundaryConditions,
                   file_name: str):
    location = boundary_conditions.location
    effective_source = boundary_conditions.effective_source
    flux_threshold = boundary_conditions.flux_threshold

    files_path = Path(Path.cwd(),
                      'simulation_files',
                      location.name,
                      effective_source.name,
                      f'Ix_{round(flux_threshold, 4)}',
                      configuration,
                      f'Nm_{number_of_mirrors}',
                      file_name)
    files_path.mkdir(parents=True, exist_ok=True)

    return files_path


def eaGeometricOptimization(number_of_mirrors: int,
                            configuration: str,
                            boundary_conditions: BoundaryConditions,
                            pop,
                            evolutionary_operators: EvolutionaryOperators,
                            file_name: str):

    # boundary conditions auxiliary variables
    bounds = boundary_conditions.bounds
    location = boundary_conditions.location
    effective_source = boundary_conditions.effective_source
    flux_threshold = boundary_conditions.flux_threshold

    # evolutionary operators auxiliary variables
    n_gens = evolutionary_operators.n_gens
    cxpb = evolutionary_operators.cxpb
    mupb = evolutionary_operators.mupb
    mustd = evolutionary_operators.mustd

    # other auxiliary variables
    cum_eff = effective_source.cum_eff
    n_dec_vars = number_of_decision_variables(number_of_mirrors=number_of_mirrors,
                                              configuration=configuration)
    #############################################################################

    # Creating the toolbox to the inputted in the GA routine
    toolbox = base.Toolbox()
    #############################################################################

    # Settings to create an individual and then a population of individuals
    toolbox.register("attr_float", random.uniform, a=0, b=1)
    toolbox.register("individual", tools.initRepeat,
                     creator.Individual, toolbox.attr_float, n=n_dec_vars)

    #############################################################################

    # Population settings
    # toolbox function for a random population
    toolbox.register("random_population",
                     tools.initRepeat,
                     list,
                     toolbox.individual)

    # toolbox function to start with a pre-defined population
    def starting_population(starting_individuals):
        # creating the starting population a list of individuals
        start_pop = [creator.Individual(ind) for ind in starting_individuals]

        # Checking if starting individuals have the proper size
        booleans = [len(ind) == n_dec_vars for ind in start_pop]
        if False in booleans:
            raise ValueError('Size of starting individuals does not match with number of decision variables')

        return start_pop

    toolbox.register("starting_population", starting_population)
    #############################################################################

    # Fitness evaluation functino of each individual
    toolbox.register("evaluate",
                     individual_evaluation,
                     configuration=configuration,
                     number_of_mirrors=number_of_mirrors,
                     bounds=bounds,
                     cum_eff=cum_eff,
                     location=location,
                     flux_threshold=flux_threshold)

    #############################################################################

    # Matting (mate) and Mutation (mutate) operators
    def mutate(ind):
        tools.mutGaussian(individual=ind, sigma=mustd, mu=0., indpb=1 / n_dec_vars)
        for i in range(len(ind)):
            ind[i] = max(0.0, min(1.0, ind[i]))

    toolbox.register('mate', tools.cxTwoPoint)
    toolbox.register('mutate', mutate)
    #############################################################################

    files_path = get_files_path(number_of_mirrors=number_of_mirrors, configuration=configuration,
                                boundary_conditions=boundary_conditions, file_name=file_name)
    evolved_population_file = Path(files_path, 'evolved_population.json')

    if evolved_population_file.is_file():
        print(f'It is {datetime.now()}! Evolved population file already exists for '
              f'{configuration} configuration, '
              f'{number_of_mirrors} mirrors, Imin = {flux_threshold} kW/m2, '
              f'{file_name}. Loading data!')
        population = get_individuals_from_file(full_file_path=evolved_population_file,
                                               return_pareto=False)
        pareto_front = get_individuals_from_file(full_file_path=evolved_population_file, return_pareto=True)
    else:
        # Population evolution routine
        print(f'It is {datetime.now()}! Starting multi-objective optimization for '
              f'{configuration} configuration, '
              f'{number_of_mirrors} mirrors, Imin = {flux_threshold} kW/m2, '
              f'{file_name}.')

        population = simpleNSGA2(pop=pop, n_gens=n_gens,
                                 cxpb=cxpb, mupb=mupb,
                                 toolbox=toolbox,
                                 files_path=files_path)

        pareto_front = [ind for ind in population if ind.rank == 1]
        hv = calculate_hypervolume(pareto_individuals=pareto_front)
        #############################################################################

        # # Exporting files
        # Data of the evolved population to be exported
        pop_dic = {'number_mirrors': number_of_mirrors,
                   'n_gens': n_gens, 'cxpb': cxpb, 'mupb0': mupb, 'mustd': mustd,
                   'configuration': configuration,
                   'population': {'individuals': population,
                                  'f1': [ind.fitness.values[0] for ind in population],
                                  'f2': [ind.fitness.values[1] for ind in population]},

                   'pareto_front': {'individuals': pareto_front,
                                    'f1': [ind.fitness.values[0] for ind in pareto_front],
                                    'f2': [ind.fitness.values[1] for ind in pareto_front],
                                    'hv': hv}
                   }

        # exporting population and
        dic2json(d=pop_dic,
                 file_path=files_path,
                 file_name=f'evolved_population')
    ################################################################################################

    return population, pareto_front


def eaSequentialOptimization(number_of_mirrors: int,
                             configuration: str,
                             boundary_conditions: BoundaryConditions,
                             eos: list,
                             gen_names: list):

    assert len(eos) == len(gen_names), 'Number of sequential generations to be evaluated does not match ' \
                                       'with number of evolutionary operators!'

    for i, (gen, eo) in enumerate(zip(gen_names, eos)):
        if i == 0:
            pop = 200
        else:
            gen_dic = read_evolved_population_file(file_name=gen_names[i - 1],
                                                   number_of_mirrors=number_of_mirrors,
                                                   configuration=configuration,
                                                   boundary_conditions=boundary_conditions)
            pop = gen_dic['population']['individuals']

        population, pareto_front = eaGeometricOptimization(number_of_mirrors=number_of_mirrors,
                                                           configuration=configuration,
                                                           boundary_conditions=boundary_conditions,
                                                           pop=pop,
                                                           evolutionary_operators=eo,
                                                           file_name=gen)

    return population, pareto_front


def get_individuals_from_file(full_file_path: Path,
                              return_pareto=False):
    if full_file_path.is_file():
        with open(full_file_path, 'r') as file:
            pop_dic = json.load(file)
    else:
        raise ValueError('This file does not exist! Please check it')

    key = 'pareto_front' if return_pareto else 'population'

    individuals = [creator.Individual(ind)
                   for ind in pop_dic[key]['individuals']]

    fitness_values = [(f1, f2)
                      for f1, f2 in zip(pop_dic[key]['f1'], pop_dic[key]['f2'])]

    for ind, fit in zip(individuals, fitness_values):
        ind.fitness.values = fit

    return individuals


def read_evolved_population_file(file_name: str,
                                 number_of_mirrors: int,
                                 configuration: str,
                                 boundary_conditions: BoundaryConditions):
    file_path = get_files_path(number_of_mirrors=number_of_mirrors, configuration=configuration,
                               boundary_conditions=boundary_conditions, file_name=file_name)
    full_file_path = Path(file_path, f'evolved_population.json')

    if full_file_path.is_file():
        with open(full_file_path, 'r') as file:
            dic = json.load(file)
    else:
        raise ValueError('This file does not exist! Please check it')

    return dic


def read_generation_file(file_name: str,
                         gen: int,
                         number_of_mirrors: int,
                         configuration: str,
                         boundary_conditions: BoundaryConditions):
    file_path = get_files_path(number_of_mirrors=number_of_mirrors, configuration=configuration,
                               boundary_conditions=boundary_conditions, file_name=file_name)
    full_file_path = Path(file_path, f'gen_{gen}.json')

    if full_file_path.is_file():
        with open(full_file_path, 'r') as file:
            dic = json.load(file)
    else:
        raise ValueError('This file does not exist! Please check it')

    return dic


def generation_hypervolume(file_name: str,
                           gen: int,
                           number_of_mirrors: int,
                           configuration: str,
                           boundary_conditions: BoundaryConditions):
    file_path = get_files_path(number_of_mirrors=number_of_mirrors, configuration=configuration,
                               boundary_conditions=boundary_conditions, file_name=file_name)
    full_file_path = Path(file_path, f'gen_{gen}.json')

    pareto_individuals = get_individuals_from_file(full_file_path=full_file_path)

    hypervolume = calculate_hypervolume(pareto_individuals=pareto_individuals)

    return hypervolume


def get_pareto_geometries(number_of_mirrors: int,
                          configuration: str,
                          boundary_conditions: BoundaryConditions,
                          file_name: str):

    print(f'It is {datetime.now()}! Importing geometries from the Pareto Front.')
    time.sleep(0.5)
    # Reading evolved population 'json' file
    population_dic = read_evolved_population_file(number_of_mirrors=number_of_mirrors,
                                                  configuration=configuration,
                                                  boundary_conditions=boundary_conditions,
                                                  file_name=file_name)

    # Getting pareto individuals and corresponding fitness values
    pareto_individuals = population_dic['pareto_front']['individuals']
    fitness_values = [(f1, f2)
                      for f1, f2 in zip(population_dic['pareto_front']['f1'],
                                        population_dic['pareto_front']['f2'])]

    # Getting the linear Fresnel geometries from the individuals "genes"
    pareto_geometries = [get_lfc_solution(individual=ind,
                                          number_of_mirrors=number_of_mirrors,
                                          configuration=configuration,
                                          bounds=boundary_conditions.bounds)

                         for ind in tqdm(pareto_individuals)]

    # Creating an attribute fitness to the lfc_geometry objects
    for geo, fit in zip(pareto_geometries, fitness_values):
        geo.fitness = fit

    return pareto_geometries


def get_pareto_geometries_data(number_of_mirrors: int,
                               configuration: str,
                               boundary_conditions: BoundaryConditions,
                               gen_name: str):

    file_path = get_files_path(number_of_mirrors=number_of_mirrors, configuration=configuration,
                               boundary_conditions=boundary_conditions, file_name=gen_name)

    full_file_path = Path(file_path, 'pareto_geometries_data.json')

    if full_file_path.is_file():
        with open(full_file_path, 'r') as opened_file:
            geo_dic = json.load(opened_file)
    else:
        # Getting the Pareto geometries with the fitness (f1, f2) of each as an attribute '.fitness'.
        pareto_geometries = get_pareto_geometries(number_of_mirrors=number_of_mirrors,
                                                  configuration=configuration,
                                                  boundary_conditions=boundary_conditions, file_name=gen_name)

        geometric_data = {
            'n': number_of_mirrors,
            'widths': [geo.primary_field.widths.tolist() for geo in pareto_geometries],
            'centers': [geo.primary_field.centers.tolist() for geo in pareto_geometries],
            'gaps': [geo.primary_field.gaps.tolist() for geo in pareto_geometries],
            'radii': [geo.primary_field.radius.tolist() for geo in pareto_geometries],
            'heights': [geo.receiver_height for geo in pareto_geometries],
            'f1': [geo.fitness[0] for geo in pareto_geometries],
            'f2': [geo.fitness[1] for geo in pareto_geometries],
            'p1': [geo.dimensionless_parameters[0] for geo in pareto_geometries],
            'p2': [geo.dimensionless_parameters[1] for geo in pareto_geometries],
            'p3': [geo.dimensionless_parameters[2] for geo in pareto_geometries]
        }

        geo_dic = dic2json(d=geometric_data,
                           file_path=file_path,
                           file_name='pareto_geometries_data')

    return geo_dic


########################################################################################################################

if __name__ == '__main__':
    pass

    # # Create the toolbox
    # toolbox = base.Toolbox()
    # # Creating the toolbox
    #
    # # Settings to create an individual and then a population of individuals
    # toolbox.register("attr_float", random.uniform, a=0, b=1)
    # toolbox.register("individual", tools.initRepeat,
    #                  creator.Individual, toolbox.attr_float, n=4)
    #
    # toolbox.register("population",
    #                  tools.initRepeat,
    #                  list,
    #                  toolbox.individual)
    #
    # ind_1 = toolbox.individual()
    # ind2 = creator.Individual([0, 0.1, 0.2, 0.3])
