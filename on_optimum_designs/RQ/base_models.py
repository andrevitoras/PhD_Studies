from on_optimum_designs.RQ.classes_and_functions import *
from on_optimum_designs.RQ.ea_functions import EvolutionaryOperators

# Optimization boundary conditions #####################################################################################
evora = SiteData(name='Evora',
                 latitude=38.53,
                 longitude=-8.0)

var_bounds = VariableBounds(height=(4 * 1e3, 20 * 1e3),
                            width=(0.2 * 1e3, 2.0 * 1e3),
                            gap=(0, 2.0 * 1e3),
                            radius=(0, 100 * 1e3))

ES1 = EffectiveSource(name='ES1',
                      sunshape='p', size=4.65e-3,
                      slope_error=2e-3, specular_error=3.e-3)

# Flux threshold condition
flux_threshold = 5.0
threshold_condition = BoundaryConditions(bounds=var_bounds,
                                         flux_threshold=flux_threshold,
                                         location=evora,
                                         effective_source=ES1)
########################################################################################################################

# Evolutionary operators
eoA = EvolutionaryOperators(n_gens=50,
                            cxpb=0.5, mupb=0.25, mustd=0.1)

eoB = EvolutionaryOperators(n_gens=50,
                            cxpb=0.8, mupb=0.4, mustd=0.2)


generations = ['genA', 'genB', 'genC', 'genD', 'genE', 'genF', 'genG', 'genH']
operators = [eoB, eoB, eoB, eoB, eoA, eoA, eoA, eoA]
####################################################################################################################

