import numpy as np
import gemmi
import rdkit

from dlib.dsmall import cif_to_mol

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
    rdkit.Chem.rdDepictor.Compute2DCoords(mol)
    coord_array = get_coord_array(mol)

    # Get bounding box
    bounds = get_bounds(mol)

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