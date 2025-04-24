from on_the_primaries_radius.base_models import *

if __name__ == '__main__':
    locations = [evora, aswan]
    geometries = [lfc_1, lfc_2]
    sources = [ES1, ES2]

    configurations = [[loc, geo, sor] for loc in locations for geo in geometries for sor in sources]

    datas = [curvature_designs_analysis(site=conf[0], geom=conf[1], source=conf[2]) for conf in configurations]
