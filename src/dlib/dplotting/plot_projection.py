import numpy as np
import gemmi
import rdkit
from rdkit.Chem import AllChem

from dlib.dsmall import cif_to_mol

def get_coord_array(mol):
    for i, conformer in enumerate(mol.GetConformers()):

        positions: np.ndarray = conformer.GetPositions()
        return positions

def get_bounds(arr, border=5.0):
    return [
        (np.min(arr, axis=0)-border)[0,1],
        (np.max(arr, axis=0)+border)[0,1]
    ]


def plot_projection(
        structure_path,
        cif_path,
        map_path
):
    # Load the structure
    st = gemmi.read_structure(str(structure_path))

    # Load the cif
    cif = gemmi.cif.read(str(cif_path))
    mol = cif_to_mol(cif)

    # Load the map
    dmap = gemmi.read_ccp4_map(str(map_path))

    # Generate the 2d projection
    AllChem.Compute2DCoords(mol)
    coord_array = get_coord_array(mol)
    print([np.min(coord_array, axis=0), np.max(coord_array, axis=0)])

    # Get bounding box
    bounds = get_bounds(mol)
    print(f"Bounds: {bounds}")

    # Generate the grid
    grid = np.meshgrid(
        np.linspace(bounds[0][0], bounds[1][0], 100),
        np.linspace(bounds[1][0], bounds[1][1], 100),

    )

    # Get the voronoi cells of the grid points relative to projection
    vcells = get_vcells(coord_array, grid)

    # For each point
        # Get the anchor atoms

        # Get the plane

        # Get the plane coords of point

        # Get the 3d plane

        # Get the 3d plane coords

        # Get the 3d pos

        # Interpolate

    # Get

    # Plot atoms

    # Plot grid points

    ...