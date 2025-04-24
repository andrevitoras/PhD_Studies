
from an_optical_method.base_models import *

files_path = Path(Path.cwd(), 'simulation_files', 'efficiency_files')
files_path.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':

    for lfc in geometries:

        [lfc.validation_analysis(ES=es, path=files_path, engine='API',
                                 force_sim=False)

         for es in effective_sources]
