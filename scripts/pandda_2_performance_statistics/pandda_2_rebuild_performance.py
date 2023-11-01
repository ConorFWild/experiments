from typing import TypedDict


import fire

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


    def get_pandda_2_rebuild_performance(self, cli=True):

        # Collect the human models

        # Collect the PanDDA 2 models

        # Get the stats for each

        ...

    def get_pandda_2_rebuild_performance_all_systems(self, cli=True):
        # Collect PanDDA 2 directories and human models

        # Loop over collecting statistics for each system

        # Combine statistics

        # Output statistics if from cli

        ...


if __name__ == "__main__":
    fire.Fire(CLI)
