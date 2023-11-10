from pathlib import Path

import fire
from scipy import stats
import pandas as pd
pd.set_option('display.max_rows', 5000)

from dlib import constants

def get_system_name_from_dtag(dtag):
    hyphens = [pos for pos, char in enumerate(dtag) if char == "-"]
    if len(hyphens) == 0:
        return None
    else:
        last_hypen_pos = hyphens[-1]
        system_name = dtag[:last_hypen_pos ]

        return system_name
def get_system_name(data_dir):
    dtags = [path.name for path in data_dir.glob('*')]
    systems = [get_system_name_from_dtag(dtag) for dtag in dtags]
    return str(stats.mode(systems, axis=None, keepdims=False)[0])

def number_of_system_hits():
    xchem_data_path = Path('/dls/labxchem/data')
    system_modelled_structures = {}
    for year_or_visit_dir in xchem_data_path.glob('*'):
        for experiment_dir in year_or_visit_dir.glob('*'):
            analysis_dir = experiment_dir / "processing" / "analysis"
            initial_models_dir = analysis_dir / "initial_models"
            model_building_dir = analysis_dir / "model_building"

            if initial_models_dir.exists():
                data_dir = initial_models_dir
            elif model_building_dir.exists():
                data_dir = model_building_dir
            else:
                print(f"No data dir for {experiment_dir}!")
                continue

            system_name = get_system_name(data_dir)

            system_modelled_structures[system_name] = 0

            for dtag_dir in data_dir.glob("*"):
                dtag = dtag_dir.name
                pandda_model_file = dtag_dir / constants.PANDDA_MODEL_FILE.format(dtag=dtag)
                if pandda_model_file.exists():
                    system_modelled_structures[system_name] += 1

    df = pd.DataFrame(
        [
            {
                "System": system_name,
                "Number of Hits": number_of_hits
            }
            for system_name, number_of_hits
            in system_modelled_structures.items()
        ]
    )
    print(df)




if __name__ == "__main__":
    fire.Fire(number_of_system_hits)