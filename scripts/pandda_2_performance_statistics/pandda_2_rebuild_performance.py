from typing import TypedDict
from pathlib import Path

import gemmi
import fire

from dlib import constants
from dlib.dcommon import read_yaml
from dlib.dsmall import match_structure_ligands, get_rmsd_from_match

class CLI:
    def compare_pandda_2_build_to_human(self,
                                        pandda_2_structure,
                                        human_structure
                                        ):
        # Get the ligand matches
        ligand_matches = match_structure_ligands(
            pandda_2_structure,
            human_structure,
        )

        # Select the closest one
        for pandda_2_lig_id, human_lig_matches in ligand_matches.items():
            distances = {}
            for human_lig_id, matched_atoms in human_lig_matches.items():
                rmsd = get_rmsd_from_match(
                    pandda_2_structure,
                    human_structure,
                    matched_atoms
                )

            closest_lig_id = min(distances, key=lambda _key: distances[_key])

        #

        ...

    def parse_human_model_dir(self, human_build_dir):
        human_structures = {}
        for path in human_build_dir.glob("*.pdb"):
            human_structures[path.stem] = gemmi.read_structure(str(path))

        return human_structures

    def parse_pandda_2_dir(self, pandda_2_dir):
        pandda_2_structures = {}
        processed_datasets_dir = pandda_2_dir / constants.PANDDA_PROCESSED_DATASETS_DIR
        for dtag_dir in processed_datasets_dir.glob('*'):
            model_building_dir = dtag_dir / constants.PANDDA_INSPECT_MODEL_DIR
            st_path = model_building_dir / constants.PANDDA_MODEL_FILE.format(dtag=dtag_dir.name)
            if st_path.exists():
                pandda_2_structures[dtag_dir.name] = gemmi.read_structure(str(st_path))

        return pandda_2_structures

    def get_pandda_2_rebuild_performance(self, pandda_2_dir, human_build_dir, cli=True):

        # Collect the human models
        print(f"\tGetting human structures...")
        human_structures = self.parse_human_model_dir(human_build_dir)
        print(human_structures)

        # Collect the PanDDA 2 models]
        print(f"\tGetting pannda 2 structures...")
        pandda_2_structures = self.parse_pandda_2_dir(pandda_2_dir)
        print(pandda_2_structures)

        # Get the stats for each

        ...

    def get_pandda_2_rebuild_performance_all_systems(self, cli=True):

        print(f"Reading datamap...")
        datamap = read_yaml("datamap.yaml")
        print(datamap)
        # Collect PanDDA 2 directories and human models

        # Loop over collecting statistics for each system
        for system, system_data in datamap.items():
            print(f"\tAnalysing system: {system}")
            self.get_pandda_2_rebuild_performance(
                Path(system_data["pandda_2"]),
                Path(system_data["human_builds"])
            )


        # Combine statistics

        # Output statistics if from cli

        ...


if __name__ == "__main__":
    fire.Fire(CLI)
