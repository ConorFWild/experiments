import os
import time

import pandas as pd

from sklearnex import patch_sklearn
patch_sklearn()

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from matplotlib import pyplot as plt
import seaborn as sns

import numpy as np
import gemmi

from rich import print as rprint

from pandda_gemmi.interfaces import *
from pandda_gemmi.args import PanDDAArgs
from pandda_gemmi.fs import PanDDAFS
from pandda_gemmi.dataset import XRayDataset, StructureArray, Structure
from pandda_gemmi.dmaps import (
    SparseDMap,
    SparseDMapStream,
    TruncateReflections,
    SmoothReflections,
)
from pandda_gemmi.alignment import Alignment, DFrame
from pandda_gemmi.processor import ProcessLocalRay, Partial
from pandda_gemmi.comparators import (
    get_comparators,
    FilterRFree,
    FilterRange,
    FilterExcludeFromAnalysis,
    FilterOnlyDatasets,
    FilterSpaceGroup,
    FilterResolution,
    FilterCompatibleStructures
)
from pandda_gemmi.event_model.event import EventBuild
from pandda_gemmi.event_model.characterization import get_characterization_sets, CharacterizationNNAndFirst
from pandda_gemmi.event_model.filter_characterization_sets import filter_characterization_sets
from pandda_gemmi.event_model.outlier import PointwiseNormal
from pandda_gemmi.event_model.cluster import ClusterDensityDBSCAN
from pandda_gemmi.event_model.score import get_model_map, ScoreCNNLigand
from pandda_gemmi.event_model.filter import (
    FilterSize,
    FilterScore,
    FilterSymmetryPosBuilds,
    FilterLocallyHighestBuildScoring
)
from pandda_gemmi.event_model.select import select_model
from pandda_gemmi.event_model.output import output_maps
from pandda_gemmi.event_model.filter_selected_events import filter_selected_events
from pandda_gemmi.site_model import HeirarchicalSiteModel, Site, get_sites
from pandda_gemmi.autobuild import AutobuildResult
from pandda_gemmi.autobuild.inbuilt import mask_dmap, get_conformers, autobuild_conformer
from pandda_gemmi.autobuild.merge import merge_autobuilds, MergeHighestBuildScore
from pandda_gemmi.ranking import rank_events, RankHighEventScore, RankHighEventScoreBySite
from pandda_gemmi.tables import output_tables
from pandda_gemmi.pandda_logging import PanDDAConsole
from pandda_gemmi import serialize

from pandda_gemmi.args import PanDDAArgs

config = {
    "Dtags": {
        "BAZ2BA-x557": {
            "Points": [
                [
                    13.11,
                    39.91,
                    29.79
                ]
            ]
        },
    },
    "Residues": [
        ("A", "1888"),
        ("A", "1914")
    ]
}
def get_sigma_map(mean_diff_array):
    # Sort array values
    sorted_mean_diff = np.sort(mean_diff_array)

    # Sample standard normal values of same size
    rng = np.random.default_rng()
    normal_samples = rng.standard_normal(sorted_mean_diff.size)

    # Sort normal values
    sorted_normal_samples = np.sort(normal_samples)

    # Get inner 50% of data
    inner_sorted_mean_diff = sorted_mean_diff[int(sorted_mean_diff.size/4):int(sorted_mean_diff.size*(3/4))]
    inner_sorted_normal_samples = sorted_normal_samples[int(sorted_mean_diff.size/4):int(sorted_mean_diff.size*(3/4))]

    # Get line of best fit angle
    z = np.polyfit(inner_sorted_mean_diff, inner_sorted_normal_samples, 1)

    return z[0]

class GetDatasetsToProcess:
    def __init__(self, filters=None):
        self.filters = filters

    def __call__(self,
                 #*args, **kwargs
                 datasets: Dict[str, DatasetInterface],
                 fs: PanDDAFSInterface
                 ):
        datasets_not_to_process = {}
        remaining_datasets = {_dtag: _dataset for _dtag, _dataset in datasets.items()}
        for _filter in self.filters:
            remaining_datasets = _filter(remaining_datasets)
            for dtag in datasets:
                if (dtag not in datasets_not_to_process) and (dtag not in remaining_datasets):
                    datasets_not_to_process[dtag] = _filter.description()

        sorted_remaining_datasets = {
            _k: remaining_datasets[_k]
            for _k
            in sorted(remaining_datasets)
        }
        sorted_datasets_not_to_process = {
            _k: datasets_not_to_process[_k]
            for _k
            in sorted(datasets_not_to_process)
        }
        return sorted_remaining_datasets, sorted_datasets_not_to_process


def main(args):
    # Record time at which PanDDA processing begins
    time_pandda_begin = time.time()

    # Create the console to print output throughout the programs run
    console = PanDDAConsole()

    # Print the PanDDA initialization message and the command line arguments
    console.start_pandda()
    console.start_parse_command_line_args()
    console.summarise_arguments(args)

    # Get the processor to handle the dispatch of functions to multiple cores and the cache of parallel
    # processed objects
    console.start_initialise_multiprocessor()
    processor: ProcessorInterface = ProcessLocalRay(args.local_cpus)
    console.print_initialized_local_processor(args)

    # Get the model of the input and output of the program on the file systems
    console.start_fs_model()
    fs: PanDDAFSInterface = PanDDAFS(Path(args.data_dirs), Path(args.out_dir), args.pdb_regex, args.mtz_regex)
    console.summarise_fs_model(fs)


    # Load the structures and reflections from the datasets found in the file system, and create references to these
    # dataset objects and the arrays of their structures in the multiprocessing cache
    console.start_load_datasets()
    datasets: Dict[str, DatasetInterface] = {
        dataset_dir.dtag: XRayDataset.from_paths(
            dataset_dir.input_pdb_file,
            dataset_dir.input_mtz_file,
            dataset_dir.input_ligands,
        )
        for dataset_dir
        in fs.input.dataset_dirs.values()
    }
    dataset_refs = {_dtag: processor.put(datasets[_dtag]) for _dtag in datasets}
    structure_array_refs = {_dtag: processor.put(StructureArray.from_structure(datasets[_dtag].structure)) for _dtag in
                            datasets}

    # Summarise the datasets loaded from the file system and serialize the information on the input into a human
    # readable yaml file
    console.summarise_datasets(datasets, fs)
    serialize.input_data(
        fs, datasets, fs.output.path / "input.yaml"
    )

    # Get the datasets to process
    datasets_to_process, datasets_not_to_process = GetDatasetsToProcess(
        [
            FilterRFree(args.max_rfree),
            FilterRange(args.dataset_range),
            FilterExcludeFromAnalysis(args.exclude_from_z_map_analysis),
            FilterOnlyDatasets(args.only_datasets)
        ]
    )(datasets, fs)
    console.summarize_datasets_to_process(datasets_to_process, datasets_not_to_process)

    # Process each dataset by identifying potential comparator datasets, constructing proposed statistical models,
    # calculating alignments of comparator datasets, locally aligning electron density, filtering statistical models
    # to the plausible set, evaluating those models for events, selecting a model to take forward based on those events
    # and outputing event maps, z maps and mean maps for that model
    pandda_events = {}
    time_begin_process_datasets = time.time()
    console.start_process_shells()
    autobuilds = {}

    # Analyse datasets
    records = []
    records_mean_diff = []
    sigma_maps = {}
    for dtag in config['Dtags']:

        # Get the dataset
        dataset = datasets[dtag]

        # Get the resolution of the dataset
        dataset_res = dataset.reflections.resolution()

        # Get the comparator datasets: these are filtered for reasonable data quality, space group compatability,
        # compatability of structural models and similar resolution
        # HACKED BUFFER TO INCLUDE ALL DATASETS OF ANY RES!
        comparator_datasets: Dict[str, DatasetInterface] = get_comparators(
            datasets,
            [
                FilterRFree(args.max_rfree),
                FilterSpaceGroup(dataset),
                FilterCompatibleStructures(dataset),
                FilterResolution(dataset_res, args.max_shell_datasets, 100, 1000)]
        )

        # Ensure the dataset itself is included in comparators
        if dtag not in comparator_datasets:
            comparator_datasets[dtag] = dataset

        # Get the resolution to process the dataset at
        processing_res = max(
            [_dataset.reflections.resolution() for _dataset in comparator_datasets.values()]
        )
        print(f"Processing resolution is: {processing_res}")

        # Print basic information about the processing to be done of the dataset
        console.begin_dataset_processing(
            dtag,
            dataset,
            dataset_res,
            comparator_datasets,
            processing_res,
            0,
            datasets_to_process,
            time_begin_process_datasets
        )

        # Get the alignments, and save them to the object store
        alignments: Dict[str, AlignmentInterface] = processor.process_dict(
            {_dtag: Partial(Alignment.from_structure_arrays).paramaterise(
                _dtag,
                structure_array_refs[_dtag],
                structure_array_refs[dtag],
            ) for _dtag in comparator_datasets}
        )
        alignment_refs = {_dtag: processor.put(alignments[_dtag]) for _dtag in comparator_datasets}

        # Get the reference frame and save it to the object store
        reference_frame: DFrame = DFrame(dataset, processor)
        reference_frame_ref = processor.put(reference_frame)

        # Get the transforms to apply to the dataset before locally aligning and save them to the object store
        transforms = [
            TruncateReflections(
                comparator_datasets,
                processing_res,
            ),
            SmoothReflections(dataset)
        ]
        transforms_ref = processor.put(transforms)

        # Load the locally aligned density maps and construct an array of them
        time_begin_get_dmaps = time.time()
        dmaps_dict = processor.process_dict(
            {
                _dtag: Partial(SparseDMapStream.parallel_load).paramaterise(
                    dataset_refs[_dtag],
                    alignment_refs[_dtag],
                    transforms_ref,
                    reference_frame_ref
                )
                for _dtag
                in comparator_datasets
            }
        )
        dmaps = np.vstack([_dmap.data.reshape((1, -1)) for _dtag, _dmap in dmaps_dict.items()])
        time_finish_get_dmaps = time.time()
        # TODO: log properly
        print(f"\t\tGot dmaps in: {round(time_finish_get_dmaps - time_begin_get_dmaps, 2)}")
        dtag_array = np.array([_dtag for _dtag in comparator_datasets])

        # Get the dataset dmap, both processed and unprocessed
        dtag_index = np.argwhere(dtag_array == dtag)
        dataset_dmap_array = dmaps[dtag_index[0][0], :]
        xmap_grid = reference_frame.unmask(SparseDMap(dataset_dmap_array))
        raw_xmap_grid = dataset.reflections.transform_f_phi_to_map(sample_rate=3)
        raw_xmap_array = np.array(raw_xmap_grid, copy=True)

        # Get the masked grid of the structure
        model_grid = get_model_map(dataset.structure.structure, xmap_grid)

        # Get the Comparator sets that define the models to try
        time_begin_get_characterization_sets = time.time()
        characterization_sets: Dict[int, Dict[str, DatasetInterface]] = get_characterization_sets(
            dtag,
            comparator_datasets,
            dmaps,
            reference_frame,
            CharacterizationNNAndFirst()
        )
        time_finish_get_characterization_sets = time.time()
        # TODO: Log properly
        print(
            f"\t\tGot characterization sets in: {round(time_finish_get_characterization_sets - time_begin_get_characterization_sets, 2)}")

        # Filter the models which are clearly poor descriptions of the density
        # In theory this step could result in the exclusion of a ground state model which provided good contrast
        # for a ligand binding in one part of the protein but fit poorly to say a large disordered region
        models_to_process, model_scores, characterization_set_masks = filter_characterization_sets(
            comparator_datasets,
            characterization_sets,
            dmaps,
            dataset_dmap_array,
            reference_frame,
            PointwiseNormal(),
        )

        # For each model, get the statistical maps
        for model_number in models_to_process:

            model_dmaps = dmaps[characterization_set_masks[model_number], :]

            # Get the statistical maps
            mean, std, z = PointwiseNormal()(
                dataset_dmap_array,
                model_dmaps
            )
            mu_grid = reference_frame.unmask(SparseDMap(mean))
            sigma_grid = reference_frame.unmask(SparseDMap(std))
            z_grid = reference_frame.unmask(SparseDMap(z))

            # Sample the distribution of the ground state datasets at the sample point
            std_array = np.array(sigma_grid)

            protein_mask_array = np.zeros(std_array.shape, dtype=bool)
            protein_mask_array[reference_frame.mask.indicies_inner_atomic] = True

            solvent_mask_array = np.zeros(std_array.shape, dtype=bool)
            solvent_mask_array[reference_frame.mask.indicies] = True
            solvent_mask_array[reference_frame.mask.indicies_inner_atomic] = False

            protein_std_array = std_array[protein_mask_array]
            model_type = f"protein_{model_number}"
            for j, val in enumerate(protein_std_array):
                records.append(
                    {
                        'Dtag': dtag,
                        'Index': j,
                        'Model': model_number,
                        'Sigma': val,
                        'ModelType': model_type
                    }
                )

            solvent_std_array = std_array[solvent_mask_array]
            model_type = f"solvent_{model_number}"
            for j, val in enumerate(solvent_std_array):
                records.append(
                    {
                        'Dtag': dtag,
                        'Index': j,
                        'Model': model_number,
                        'Sigma': val,
                        'ModelType': model_type
                    }
                )


            #####################

            mean_difference_array = np.array(xmap_grid) - np.array(mu_grid)

            protein_mean_diff_array = mean_difference_array[protein_mask_array]
            model_type = f"mean_diff_protein_{model_number}"
            for j, val in enumerate(protein_mean_diff_array):
                records_mean_diff.append(
                    {
                        'Dtag': dtag,
                        'Index': j,
                        'Model': model_number,
                        'MeanDiff': val,
                        'ModelType': model_type
                    }
                )
            sigma_map = get_sigma_map(protein_mean_diff_array)
            sigma_maps[model_type] = sigma_map

            solvent_mean_diff_array = mean_difference_array[solvent_mask_array]
            model_type = f"mean_diff_solvent_{model_number}"
            for j, val in enumerate(solvent_mean_diff_array):
                records_mean_diff.append(
                    {
                        'Dtag': dtag,
                        'Index': j,
                        'Model': model_number,
                        'MeanDiff': val,
                        'ModelType': model_type
                    }
                )
            sigma_map = get_sigma_map(solvent_mean_diff_array)
            sigma_maps[model_type] = sigma_map



        # Plot in seaborn
        data = pd.DataFrame(records)
        fig, ax = plt.subplots(
            figsize=(6, 4.8)
        )

        sns.ecdfplot(
            data=data,
            ax=ax,
            x='Sigma',
            hue='ModelType'
        )
        ax.set_xscale('log')

        # Save
        plt.savefig('outputs/std_protein_solvent_dist.png')

        data.to_csv('outputs/std_protein_solvent_dist.csv')


        # Plot in seaborn
        data = pd.DataFrame(records_mean_diff)
        fig, ax = plt.subplots(
            figsize=(6, 4.8)
        )

        sns.ecdfplot(
            data=data,
            ax=ax,
            x='MeanDiff',
            hue='ModelType'
        )
        # ax.set_xscale('log')

        # Save
        plt.savefig('outputs/mean_diff_protein_solvent_dist.png')

        data.to_csv('outputs/mean_diff_protein_solvent_dist.csv')

        rprint(sigma_maps)

if __name__ == '__main__':
    # Parse Command Line Arguments
    args = PanDDAArgs.from_command_line()

    # Process the PanDDA
    main(args)
