import numpy as np
import gemmi
from scipy.spatial import KDTree
import rdkit
from rdkit.Chem import AllChem

from dlib.dsmall import cif_to_mol

def get_coord_array(mol):
    for i, conformer in enumerate(mol.GetConformers()):

        positions: np.ndarray = conformer.GetPositions()
        return positions[:, np.array([0,1])]

def get_bounds(arr, border=5.0):
    return [
        (np.min(arr, axis=0)-border)[np.array([0,1]),],
        (np.max(arr, axis=0)+border)[np.array([0,1]),]
    ]

def get_vcells(kd, grid_samples):





    # Get nearest neighbour
    nbs = kd.query(grid_samples)

    return nbs[1]

def get_structure_atom_array(st):
    atoms = {}
    for model in st:
        for chain in model:
            for res in chain:
                if res.name == "LIG":
                    for atom in res:
                        pos = atom.pos
                        atoms[atom.name] = [
                            pos.x,
                            pos.y,
                            pos.z
                        ]

    return atoms


def plot_projection(
        structure_path,
        cif_path,
        map_path
):
    # Load the structure
    st = gemmi.read_structure(str(structure_path))

    # Structure atom array
    st_atom_pos_dict = get_structure_atom_array(st)

    # Load the cif
    cif = gemmi.cif.read(str(cif_path))
    mol, atom_ids = cif_to_mol(cif)

    # # Get the atom id array
    # atom_id_array = get_atom_id_array(mol)

    # Load the map
    dmap = gemmi.read_ccp4_map(str(map_path))

    # Generate the 2d projection
    AllChem.Compute2DCoords(mol)
    coord_array = get_coord_array(mol)
    print([np.min(coord_array, axis=0), np.max(coord_array, axis=0)])

    # Get bounding box
    bounds = get_bounds(coord_array)
    print(f"Bounds: {bounds}")

    # Generate the grid
    grid = np.meshgrid(
        np.linspace(bounds[0][0], bounds[1][0], 100),
        np.linspace(bounds[1][0], bounds[1][1], 100),
    )
    # Get the query array
    xs = grid[0].flatten()
    ys = grid[1].flatten()
    grid_samples = np.hstack(
        [
            xs.reshape(-1,1),
            ys.reshape(-1, 1),
        ]
    )
    # Get the voronoi cells of the grid points relative to projection
    # Get the atom coord kdtree
    kd = KDTree(
        coord_array
    )
    vcells = get_vcells(kd, grid_samples)
    print(vcells)

    # For each point
    for sample in grid_samples:
    # for j, atom_id in enumerate(atom_ids):
        print(sample)
        # Get the anchor atoms
        nbs = kd.query(sample, k=3)[1]
        print(nbs)
        print(np.array(atom_ids)[nbs])
        # exit()

        # Get structure poss
        nbr_poss = {}
        for nbr in nbs:
            pos = st_atom_pos_dict[atom_ids[nbr]]
            nbr_poss[atom_ids[nbr]] = pos
        print(nbr_poss)
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