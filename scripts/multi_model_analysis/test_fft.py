import numpy as np
import gemmi
import itertools
from scipy.spatial.transform import Rotation as R
from scipy.signal import oaconvolve
import time


def get_xmap(mtz_path, structure_array, step=0.5, border=6.0):
    mtz = gemmi.read_mtz_file(mtz_path)

    grid = mtz.transform_f_phi_to_map('FWT', 'PHWT', sample_rate=3.0)

    # Get homogenous grid
    mini, maxi = np.min(structure_array, axis=0) - border, np.max(structure_array, axis=0) + border
    diff = maxi - mini
    gridding = np.ceil(diff / step)
    print(gridding)

    # then we setup a transformation (array indices) -> (position [A]).

    tr = gemmi.Transform()
    tr.mat.fromlist([[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]])
    tr.vec.fromlist([mini[0], mini[1], mini[2]])

    # finally we calculate interpolated values
    arr = np.zeros((int(gridding[0]), int(gridding[1]), int(gridding[2])), dtype=np.float32)
    grid.interpolate_values(arr, tr)
    new_grid = gemmi.FloatGrid(int(gridding[0]), int(gridding[1]), int(gridding[2]))
    new_grid.set_unit_cell(gemmi.UnitCell(gridding[0] * 0.5, gridding[1] * 0.5, gridding[2] * 0.5, 90.0, 90.0, 90.0))
    new_grid_array = np.array(new_grid, copy=False)
    new_grid_array[:, :, :] = arr[:, :, :]

    return new_grid

    ...

def get_structure_array(structure_path):
    st = gemmi.read_structure(structure_path)
    poss = []
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    pos = atom.pos
                    poss.append([pos.x,pos.y,pos.z])

    return np.array(poss)


def get_structure_map(structure_array, x, y, z, step=0.5, border=2.0):
    centroid = np.mean(structure_array, axis=0)
    r = R.from_euler('xyz', [x, y, z], degrees=True)
    centered_array = structure_array - centroid

    rotated_array = r.apply(centered_array)
    mini, maxi = np.min(rotated_array, axis=0) - border, np.max(rotated_array, axis=0) + border
    diff = maxi - mini
    new_centroid = diff / 2
    new_array = rotated_array + new_centroid

    print(np.min(new_array, axis=0))

    gridding = np.ceil(diff / step)
    new_grid = gemmi.FloatGrid(
        int(gridding[0]),
        int(gridding[1]),
        int(gridding[2])
    )
    new_grid.set_unit_cell(gemmi.UnitCell(gridding[0] * 0.5, gridding[1] * 0.5, gridding[2] * 0.5, 90.0, 90.0, 90.0))

    for row in new_array:
        pos = gemmi.Position(row[0], row[1], row[2])
        new_grid.set_points_around(pos, 1.0, 1.0)

    return new_grid

    ...

def fft_convolve(ligand_map, xmap):
    ligand_map_array = np.array(ligand_map, copy=False)
    xmap_array = np.array(xmap, copy=False)
    conv = oaconvolve(xmap_array, ligand_map_array, mode='same')
    #print([ligand_map_array.shape, xmap_array.shape, conv.shape])
    return conv
    ...


# Define the structure and xmap paths
mtz_path = '/dls/labxchem/data/2020/lb18145-153/processing/analysis/model_building/Mpro-x0089/dimple.mtz'
protein_structure_path = '/dls/labxchem/data/2020/lb18145-153/processing/analysis/model_building/Mpro-x0089/dimple.pdb'
ligand_structure_path = '/dls/labxchem/data/2020/lb18145-153/processing/analysis/model_building/Mpro-x0089/Z1650168321.pdb'

# Get the structure map
protein_structure_array = get_structure_array(protein_structure_path)
ligand_structure_array = get_structure_array(ligand_structure_path)

# Get aligned structure
mini = np.min(protein_structure_array, axis=0)
st = gemmi.read_structure(protein_structure_path)
for model in st:
    for chain in model:
        for res in chain:
            for atom in res:
                pos = atom.pos
                pos.x = pos.x - (mini[0] - 6.0)
                pos.y = pos.y - (mini[1] - 6.0)
                pos.z = pos.z - (mini[2] - 6.0)
st.write_minimal_pdb('out.pdb')

# Get the xmap
xmap = get_xmap(mtz_path, protein_structure_array)
ccp4 = gemmi.Ccp4Map()
ccp4.grid = xmap
ccp4.update_ccp4_header()
ccp4.write_ccp4_map('out.ccp4')

for x, y, z in itertools.product(np.linspace(0, 360, 10, endpoint=False), np.linspace(0, 360, 10, endpoint=False),
                                 np.linspace(0, 360, 10, endpoint=False)):
    begin_fft = time.time()
    structure_map = get_structure_map(structure_array, x, y, z)
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = structure_map
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map('out_ligand.ccp4')
    fft = fft_convolve(structure_map, xmap)
    finish_fft = time.time()
    print(f'FFTd in {finish_fft - begin_fft}')

    new_grid = gemmi.FloatGrid(
        fft.shape[0],
        fft.shape[1],
        fft.shape[2]
    )
    new_grid.set_unit_cell(gemmi.UnitCell(fft.shape[0] * 0.5, fft.shape[1] * 0.5, fft.shape[2] * 0.5, 90.0, 90.0, 90.0))
    new_grid_array = np.array(new_grid, copy=False)
    new_grid_array[:, :, :] = fft[:, :, :]
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = new_grid
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map('out_fft.ccp4')

    raise Exception

