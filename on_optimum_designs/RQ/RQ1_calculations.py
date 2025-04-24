
from on_optimum_designs.RQ.base_models import *
from on_optimum_designs.RQ.ea_functions import *


###################################################################################


if __name__ == '__main__':

    mirrors_numbers = [12, 25, 30]
    for n in mirrors_numbers:

        _, _ = eaSequentialOptimization(
            number_of_mirrors=n,
            configuration='non-uniform',
            boundary_conditions=threshold_condition,
            gen_names=generations,
            eos=operators)

        _, _ = eaSequentialOptimization(
            number_of_mirrors=n,
            configuration='variable-radius',
            boundary_conditions=threshold_condition,
            gen_names=generations,
            eos=operators)
