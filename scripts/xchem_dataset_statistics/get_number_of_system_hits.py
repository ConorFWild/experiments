from pathlib import Path

import fire
import numpy as np
from scipy import stats
import pandas as pd
pd.set_option('display.max_rows', 5000)
pd.set_option('display.max_colwidth', 100)
pd.set_option('display.width', 200)


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
    if len(dtags) < 2:
        return None
    # print(f"Got {len(dtags)} datasets")
    systems = [get_system_name_from_dtag(dtag) for dtag in dtags ]
    systems_sanitized = [system for system in systems if system]
    if len(systems_sanitized) < 2:
        return None
    unique_vals, counts = np.unique(systems_sanitized, return_counts=True)

    mode_idx = np.argmax(counts)
    mode = unique_vals[mode_idx]

    return str(mode)

def number_of_system_hits():
    xchem_data_path = Path('/dls/labxchem/data')
    records = []
    for year_or_visit_dir in xchem_data_path.glob('*'):
        for experiment_dir in year_or_visit_dir.glob('*'):
            print(f"Processing: {experiment_dir}")
            analysis_dir = experiment_dir / "processing" / "analysis"
            initial_models_dir = analysis_dir / "initial_model"
            model_building_dir = analysis_dir / "model_building"

            if initial_models_dir.exists():
                data_dir = initial_models_dir
            elif model_building_dir.exists():
                data_dir = model_building_dir
            else:
                print(f"No data dir for {experiment_dir}!")
                continue

            try:
                print(f"Got data dir: {data_dir}")
                system_name = get_system_name(data_dir)
                if system_name is None:
                    print(f"No datasets! Skipping!")
                    continue

                num_modelled_structures = 0

                for dtag_dir in data_dir.glob("*"):
                    dtag = dtag_dir.name
                    pandda_model_file = dtag_dir / constants.PANDDA_MODEL_FILE.format(dtag=dtag)
                    if pandda_model_file.exists():
                        num_modelled_structures += 1


                records.append(
                    {
                        "System": system_name,
                        "Number of Hits": num_modelled_structures,
                        "Data Dir": data_dir
                    }
                )

            except Exception as e:
                print(f"An exception occured!")
                print(e)

    df = pd.DataFrame(
        records
    )
    print(df)




if __name__ == "__main__":
    fire.Fire(number_of_system_hits)