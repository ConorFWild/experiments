from typing import TypedDict
from pathlib import Path

import gemmi
import fire
import pandas as pd

from dlib import constants
from dlib.dcommon import read_yaml
from dlib.dsmall import match_structure_ligands, get_rmsd_from_match, get_rmsd_from_closest_atom

class CLI:
    def compare_pandda_2_build_to_human(self,
                                        pandda_2_structure,
                                        human_structure
                                        ):
        # Get the ligand matches
        print(f"Getting ligand matches...")
        ligand_matches = match_structure_ligands(
            pandda_2_structure,
            human_structure,
        )
        # print(f"Ligand matches: {ligand_matches}")

        # Select the closest one
        distances_atom_match = {}
        for match_id, matched_atoms in ligand_matches.items():
            # distances = {}
            # for , matched_atoms in human_lig_matches.items():
            rmsd = get_rmsd_from_match(
                pandda_2_structure,
                human_structure,
                matched_atoms
            )
            # print(f"\t\tRMSD is: {rmsd}")
            distances_atom_match[match_id] = rmsd

        closest_atom_match_match_id = min(
            distances_atom_match,
            key=lambda _key: distances_atom_match[_key],
        )
        # print(f"Closest match id is: {closest_match_id} : {distances[closest_match_id]}")

        # Get the RMSD by closest atom
        rmsd_closest_atom = get_rmsd_from_closest_atom(
            pandda_2_structure,
            human_structure,
            closest_atom_match_match_id[0],
            closest_atom_match_match_id[1]
        )

        return {
            "PanDDA 2 LIG ID": str(closest_atom_match_match_id[0]),
            "Human Build LIG ID": str(closest_atom_match_match_id[1]),
            "RMSD (atom match)": distances_atom_match[closest_atom_match_match_id],
            "RMSD (closest atom)":
        }

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
        # print(human_structures)

        # Collect the PanDDA 2 models]
        print(f"\tGetting pannda 2 structures...")
        pandda_2_structures = self.parse_pandda_2_dir(pandda_2_dir)
        # print(pandda_2_structures)

        # Get the stats for each
        stats = {}
        for dtag, human_build in human_structures.items():
            if dtag not in pandda_2_structures:
                dtag_stats = {

                        "PanDDA 2 LIG ID": None,
                        "Human Build LIG ID": None,
                        "RMSD": None

                }
            else:
                dtag_stats = self.compare_pandda_2_build_to_human(
                    pandda_2_structures[dtag],
                    human_build
                )
            stats[dtag] = dtag_stats

        return stats

        ...

    def get_pandda_2_rebuild_performance_all_systems(self, cli=True):

        print(f"Reading datamap...")
        datamap = read_yaml("datamap.yaml")
        print(datamap)
        # Collect PanDDA 2 directories and human models

        # Loop over collecting statistics for each system
        stats = {}
        for system, system_data in datamap.items():
            print(f"\tAnalysing system: {system}")
            system_stats = self.get_pandda_2_rebuild_performance(
                Path(system_data["pandda_2"]),
                Path(system_data["human_builds"])
            )
            stats[system] = system_stats


        # Combine statistics
        table = pd.DataFrame(
            [
                {
                    "System": system,
                    "Dtag": dtag,
                    "PanDDA 2 LIG ID": dtag_stats['PanDDA 2 LIG ID'],
            "Human Build LIG ID": dtag_stats['Human Build LIG ID'],
            "RMSD": dtag_stats['RMSD']
                }
                for system, system_stats in stats.items()
                for dtag, dtag_stats in system_stats.items()
            ]
        )
        print(table)

        # Output statistics if from cli

        ...


if __name__ == "__main__":
    fire.Fire(CLI)
