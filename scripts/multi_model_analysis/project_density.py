from pathlib import Path

import fire

from dlib.dplotting.plot_projection import plot_projection


data = {
    "BAZ2BA":
        {
            "ModelDir": Path('../../Downloads/pandda_BAZ2BA_all_hits_all_models/processed_datasets'),
            "Datasets":
                {
                    "BAZ2BA-x446": {
                        "ModelPath": "autobuild/1_4_ligand_0.pdb",
                        "CifPath": 'ligand_files/ligand.cif',
                    },
                    "BAZ2BA-x557": {
                        "ModelPath": "autobuild/2_4_ligand_0.pdb",
                        "CifPath": 'ligand_files/ligand.cif'
                    },
                    "BAZ2BA-x529": {
                        "ModelPath": "autobuild/3_7_ligand_0.pdb",
                        "CifPath": 'ligand_files/ligand.cif'
                    },

                }
        }
}

def main(
        # model_path,
        # cif_path,
        # map_path
):
    #

    base_path = Path('outputs')

    for system, system_data in data.items():
        print(system)
        model_dir = system_data['ModelDir']
        for dtag, dtag_data in system_data['Datasets'].items():
            print(f"\t{dtag}")

            #

            dtag_dataset = model_dir / dtag
            model_path = dtag_dataset / dtag_data['ModelPath']
            cif_path = dtag_dataset /dtag_data['CifPath']

            for map_path in (dtag_dataset / 'model_maps').glob('*'):
                print(f'\t\t{map_path.name}')

                # Produce and save figure for default mu map
                output_path = base_path / f'{dtag}_{map_path.stem}.png'
                plot_projection(
                    model_path,
                    cif_path,
                    map_path,
                    output_path
                )

            # Produce and save figure for pandda2 mu map

            # Produce and save figure for default sigma map

            # Produce and save figure for pandda2 sigma map

            # Produce and save figure for default z map

            # Produce and save figure for pandda2 z map


    ...

if __name__ == "__main__":
    fire.Fire(main)