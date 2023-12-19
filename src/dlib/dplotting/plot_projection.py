import itertools

import numpy as np
import gemmi
from scipy.spatial import KDTree
from scipy import spatial
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, RBFInterpolator
import rdkit
from rdkit.Chem import AllChem
from matplotlib import pyplot as plt
from sklearn import decomposition

from dlib.dsmall import cif_to_mol

from rich import print as rprint

def get_coord_array(mol):
    for i, conformer in enumerate(mol.GetConformers()):

        positions: np.ndarray = conformer.GetPositions()
        return positions[:, np.array([0,1])]

def get_bounds(arr, border=3.0):
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
                        atoms[atom.name] = np.array([
                            pos.x,
                            pos.y,
                            pos.z
                        ])

    return atoms

def get_plane_vectors(nbr_poss):
    pvs = {}
    atom_ids = list(nbr_poss.keys())

    pv1 = nbr_poss[atom_ids[0]] - nbr_poss[atom_ids[1]]
    pvs[(atom_ids[1], atom_ids[0])] = pv1 #/ np.linalg.norm(pv1)

    pv2 = nbr_poss[atom_ids[2]] - nbr_poss[atom_ids[1]]
    pvs[(atom_ids[1], atom_ids[2])] = pv2 #/ np.linalg.norm(pv2)

    return pvs

# def plot_projection(
#         structure_path,
#         cif_path,
#         map_path
# ):
#     # Load the structure
#     st = gemmi.read_structure(str(structure_path))
#
#     # Structure atom array
#     st_atom_pos_dict = get_structure_atom_array(st)
#
#     # Load the cif
#     cif = gemmi.cif.read(str(cif_path))
#     mol, atom_ids = cif_to_mol(cif)
#
#     # # Get the atom id array
#     # atom_id_array = get_atom_id_array(mol)
#
#     # Load the map
#     ccp4 = gemmi.read_ccp4_map(str(map_path), )
#     ccp4.setup(0.0)
#     dmap = ccp4.grid
#
#     # Generate the 2d projection
#     AllChem.Compute2DCoords(mol)
#     coord_array = get_coord_array(mol)
#     print(coord_array)
#     print([np.min(coord_array, axis=0), np.max(coord_array, axis=0)])
#
#     # Get bounding box
#     bounds = get_bounds(coord_array)
#     print(f"Bounds: {bounds}")
#
#     # Generate the grid
#     # grid = np.meshgrid(
#     #     np.linspace(bounds[0][0], bounds[1][0], 100),
#     #     np.linspace(bounds[1][0], bounds[1][1], 100),
#     # )
#     # Get the query array
#     # xs = grid[0].flatten()
#     # ys = grid[1].flatten()
#     xs = np.linspace(bounds[0][0], bounds[1][0], 500)
#     ys = np.linspace(bounds[0][1], bounds[1][1], 500)
#     grid_samples = np.array(
#         [x for x in itertools.product(
#             xs, ys
#         )]
#     )
#     print(grid_samples)
#     # exit()
#     print(grid_samples.shape)
#     # Get the voronoi cells of the grid points relative to projection
#     # Get the atom coord kdtree
#     kd = KDTree(
#         coord_array
#     )
#     vcells = get_vcells(kd, grid_samples)
#     print(vcells)
#
#     # exit()
#
#     # For each point
#     values = []
#     for sample in grid_samples:
#     # for j, atom_id in enumerate(atom_ids):
#     #     print(sample)
#         # Get the anchor atoms
#         dists, nbs = kd.query(sample, k=3)
#         if dists[0] > 3:
#             # print(dists[0])
#             values.append(0)
#             continue
#         # print(nbs)
#         # print(np.array(atom_ids)[nbs])
#         # exit()
#
#         # Get the 2d Poss
#         nbr_poss = {}
#         for nbr in nbs:
#             pos = coord_array[nbr]
#             nbr_poss[atom_ids[nbr]] = pos
#         # print(nbr_poss)
#         # Get the plane
#         pvs = get_plane_vectors(nbr_poss)
#         # print(pvs)
#
#         # Get the plane coords of point relative to p1
#         pv_keys = list(pvs.keys())
#         point_rel = sample - nbr_poss[pv_keys[0][0]]
#         pv1 = pvs[pv_keys[0]]
#         # comp1 =  np.dot(point_rel, pv1) / np.linalg.norm(pv1)
#
#         pv2 = pvs[pv_keys[1]]
#         # comp2 =  np.dot(point_rel, pv2) / np.linalg.norm(pv2)
#         mat = np.vstack(
#                 [
#                     pv1.reshape(1,2),
#                     pv2.reshape(1,2)
#                 ])
#         components = np.linalg.solve(
#             mat.T,
#             point_rel
#         )
#
#
#
#         # Get structure poss
#         nbr_poss_3d = {}
#         for nbr in nbs:
#             pos = st_atom_pos_dict[atom_ids[nbr]]
#             nbr_poss_3d[atom_ids[nbr]] = pos
#
#         # Get the 3d plane
#         pvs_3d = get_plane_vectors(nbr_poss_3d)
#         pv1_3d = pvs_3d[pv_keys[0]]
#         pv2_3d = pvs_3d[pv_keys[1]]
#
#         # Get the 3d pos
#         point_3d_rel = (components[0] * pv1_3d) + (components[1] * pv2_3d)
#         point_3d = point_3d_rel + nbr_poss_3d[pv_keys[0][0]]
#
#         # Interpolate
#         value = dmap.interpolate_value(
#         # value=dmap.tricubic_interpolation(
#                 gemmi.Position(
#                 point_3d[0],
#                 point_3d[1],
#                 point_3d[2],
#             )
#         )
#         values.append(
#             value
#         )
#
#
#         # if np.linalg.norm(point_3d_rel) < 0.5:
#         #     rprint({
#         #         "Pos": sample,
#         #         "Anchor Poss": nbr_poss,
#         #         "Relative Pos": point_rel,
#         #         "Relative Pos Distance": np.linalg.norm(point_rel),
#         #         "Plane Vector 1": pv1,
#         #         "Plane Vector 2": pv2,
#         #         "Components": components,
#         #         "Reconstruction": (components[0] * pv1) + (components[1] * pv2),
#         #         "3D Plane Vectors": pvs_3d,
#         #         "Point 3D": point_3d,
#         #         "Anchor Poss 3D": nbr_poss_3d,
#         #         "Point 3d Relative": point_3d_rel,
#         #         "Point 3d Relative Dist": np.linalg.norm(point_3d_rel),
#         #         # "Reconstruction 2": np.dot(mat, components)
#         #         "Value": value
#         #     })
#             # exit()
#
#
#     # h = plt.contourf(xs, ys, np.array(values).reshape(100,100))
#     # h = plt.imshow(np.array(values).reshape(100,100))
#     plt.figure(figsize=(16, 9))
#
#     h = plt.scatter(
#         grid_samples[:,0],
#         grid_samples[:,1],
#         c=values
#     )
#     plt.scatter(
#         coord_array[:,0],
#         coord_array[:,1],
#         c='r'
#     )
#     plt.axis('scaled')
#     # plt.colorbar()
#     print(f"Writing map!")
#     plt.savefig('test.png')
#
#
#         # exit()
#
#     # Get
#
#     # Plot atoms
#
#     # Plot grid points
#
#     ...
#
# def get_transform(
#         ref,
#         mov
# ):
#
#         mean_mov = np.mean(mov, axis=0)
#         mean_ref = np.mean(ref, axis=0)
#
#         vec = np.array([0.0, 0.0, 0.0])
#
#         de_meaned_mov = mov - mean_mov
#         de_meaned_ref = ref - mean_ref
#
#         rotation, rmsd = spatial.transform.Rotation.align_vectors(de_meaned_mov, de_meaned_ref)
#
#         return vec, rotation.as_euler('zyx'), mean_ref, mean_mov
#
# def plot_projection(structure_path,
#         cif_path,
#         map_path
# ):
#     # Load the structure
#     st = gemmi.read_structure(str(structure_path))
#
#     # Structure atom array
#     st_atom_pos_dict = get_structure_atom_array(st)
#
#     # Load the cif
#     cif = gemmi.cif.read(str(cif_path))
#     mol, atom_ids, atom_type_loop = cif_to_mol(cif)
#
#     # # Get the atom id array
#     # atom_id_array = get_atom_id_array(mol)
#
#     # Load the map
#     ccp4 = gemmi.read_ccp4_map(str(map_path), )
#     ccp4.setup(0.0)
#     dmap = ccp4.grid
#
#     # Generate the 2d projection
#     AllChem.Compute2DCoords(mol)
#     mask = []
#     for el in atom_type_loop:
#         if el == "H":
#             mask.append(False)
#         else:
#             mask.append(True)
#     coord_array = get_coord_array(mol)[np.array(mask), :]
#     print(coord_array)
#     print([np.min(coord_array, axis=0), np.max(coord_array, axis=0)])
#
#     # Get bounding box
#     bounds = get_bounds(coord_array)
#     print(f"Bounds: {bounds}")
#
#     # Generate the grid
#
#     xs = np.linspace(bounds[0][0], bounds[1][0], 500)
#     ys = np.linspace(bounds[0][1], bounds[1][1], 500)
#     grid_samples = np.array(
#         [x for x in itertools.product(
#             xs, ys
#         )]
#     )
#     print(grid_samples)
#     print(grid_samples.shape)
#     # Get the voronoi cells of the grid points relative to projection
#     # Get the atom coord kdtree
#     kd = KDTree(
#         coord_array
#     )
#
#     # For each atom get the three closest
#     transforms = {}
#     for j, pos in enumerate(coord_array):
#         dists, nbs = kd.query(pos, k=3)
#         poss_2d_unpadded = np.array(
#             [
#                 coord_array[nbr]
#                 for nbr
#                 in nbs
#                 ]
#             )
#
#         print(poss_2d_unpadded)
#
#         poss_2d = np.pad(
#             poss_2d_unpadded,
#             [(0, 0), (0, 1)]
#         )
#         poss_3d = np.array(
#             [
#             st_atom_pos_dict[atom_ids[nbr]]
#             for nbr
#             in nbs
#             ]
#         )
#         print(poss_2d)
#         print(poss_3d)
#         transform = get_transform(
#             poss_2d,
#             poss_3d
#         )
#         rprint(transform)
#         transforms[atom_ids[j]] = transform
#     for tr in transforms.values():
#         print(tr[1])
#
#     interpoland = np.array(
#         [
#             np.concatenate(
#                 [
#                     # tr[1].flatten(),
#                     # tr[2].flatten(),
#                     tr[3].flatten()
#                 ]
#             )
#             for tr
#             in transforms.values()
#             ]
#         )
#     print(interpoland)
#
#     # transform_interpolator = LinearNDInterpolator(
#     #     coord_array,
#     #     interpoland
#     # )
#     transform_interpolator = RBFInterpolator(
#         coord_array,
#         interpoland,
#         smoothing=2
#     )
#     # transform_interpolator_nearest = NearestNDInterpolator(
#     #     coord_array,
#     #     interpoland
#     # )
#
#     # For each point
#     values = []
#     for sample in grid_samples:
#         # for j, atom_id in enumerate(atom_ids):
#         #     print(sample)
#         # Get the anchor atoms
#         dists, nbs = kd.query(sample, k=3)
#         if dists[0] > 2:
#             # print(dists[0])
#             values.append(0)
#             continue
#         # print(nbs)
#         # print(np.array(atom_ids)[nbs])
#         # exit()
#
#         # Get the 2d Poss
#         nbr_poss = {}
#         for nbr in nbs:
#             pos = coord_array[nbr]
#             nbr_poss[atom_ids[nbr]] = pos
#
#         # Get the transform
#         # tr = transforms[atom_ids[nbs[0]]]
#         # tr_array = transform_interpolator(sample)[0]
#         tr_array = transform_interpolator(sample.reshape(1,2))[0]
#
#         # print(f"Transform array")
#         # print(tr_array)
#         # if np.isnan(tr_array[0]):
#         # tr_array = transform_interpolator_nearest(sample)[0]
#             # print(tr_array)
#
#         # mat = tr_array[:9].reshape(3,3)
#         # ref = tr_array[9:12]
#         # mov = tr_array[12:16]
#         # rot = tr_array[:3]
#         # ref = tr_array[3:6]
#         # mov = tr_array[6:9]
#
#         # print(mat)
#         # print(ref)
#         # print(mov)
#
#
#         # Get the sample point
#         # sample_point_2d = np.array([sample[0], sample[1], 0.0])
#         # sample_point_2d_rel = sample_point_2d - ref
#         # sample_point_3d_rel = np.matmul(
#         #     spatial.transform.Rotation.from_euler('zyx', rot).as_matrix(),
#         #     sample_point_2d_rel,
#         # )
#         # point_3d = sample_point_3d_rel + mov
#
#         point_3d = tr_array
#
#         # Interpolate
#         value = dmap.interpolate_value(
#             # value=dmap.tricubic_interpolation(
#             gemmi.Position(
#                 point_3d[0],
#                 point_3d[1],
#                 point_3d[2],
#             )
#         )
#         values.append(
#             value
#         )
#
#         if dists[0] < 1.0:
#             rprint({
#             "Pos": sample,
#             # "Anchor Pos 2d": tr[2],
#             # "Anchor Pos 3d": tr[3],
#             # "Relative Pos": sample_point_2d_rel,
#             # "Relative Pos Distance": np.linalg.norm(sample_point_2d_rel),
#             "Point 3D": point_3d,
#             # "Point 3d Relative": sample_point_3d_rel,
#             # "Point 3d Relative Dist": np.linalg.norm(sample_point_3d_rel),
#             # "Reconstruction 2": np.dot(mat, components)
#             "Value": value
#         })
#
#
#     plt.figure(figsize=(16, 9))
#
#     h = plt.scatter(
#         grid_samples[:,0],
#         grid_samples[:,1],
#         c=values
#     )
#     plt.scatter(
#         coord_array[:,0],
#         coord_array[:,1],
#         c='r'
#     )
#     plt.axis('scaled')
#     plt.colorbar()
#     print(f"Writing map!")
#     plt.savefig('test.png')
#
# # def plot_projection(structure_path,
# #         cif_path,
# #         map_path
# # ):
# #     # Load the structure
# #     st = gemmi.read_structure(str(structure_path))
# #
# #     # Structure atom array
# #     st_atom_pos_dict = get_structure_atom_array(st)
# #
# #     # Load the cif
# #     cif = gemmi.cif.read(str(cif_path))
# #     mol, atom_ids, atom_type_loop = cif_to_mol(cif)
# #
# #     # # Get the atom id array
# #     # atom_id_array = get_atom_id_array(mol)
# #
# #     # Load the map
# #     ccp4 = gemmi.read_ccp4_map(str(map_path), )
# #     ccp4.setup(0.0)
# #     dmap = ccp4.grid
# #
# #     # Generate the 2d projection
# #     AllChem.Compute2DCoords(mol)
# #     mask = []
# #     for el in atom_type_loop:
# #         if el == "H":
# #             mask.append(False)
# #         else:
# #             mask.append(True)
# #     coord_array = get_coord_array(mol)[np.array(mask), :]
# #     print(coord_array)
# #     print([np.min(coord_array, axis=0), np.max(coord_array, axis=0)])
# #
# #     # Get bounding box
# #     bounds = get_bounds(coord_array)
# #     print(f"Bounds: {bounds}")
# #
# #     # Generate the grid
# #
# #     xs = np.linspace(bounds[0][0], bounds[1][0], 500)
# #     ys = np.linspace(bounds[0][1], bounds[1][1], 500)
# #     grid_samples = np.array(
# #         [x for x in itertools.product(
# #             xs, ys
# #         )]
# #     )
# #     print(grid_samples)
# #     print(grid_samples.shape)
# #     # Get the voronoi cells of the grid points relative to projection
# #     # Get the atom coord kdtree
# #     kd = KDTree(
# #         coord_array
# #     )
# #
# #     # For each atom get the three closest
# #     planes = {}
# #     for j, pos in enumerate(coord_array):
# #         dists, nbs = kd.query(pos, k=3)
# #
# #         # Get the 2d poss
# #         poss_2d = np.array(
# #             [
# #                 coord_array[nbr]
# #                 for nbr
# #                 in nbs
# #                 ]
# #             )
# #
# #         # Get the 3d poss
# #         poss_3d = np.array(
# #             [
# #             st_atom_pos_dict[atom_ids[nbr]]
# #             for nbr
# #             in nbs
# #             ]
# #         )
# #
# #         # Get the 2d plane vecs centered on atom
# #
# #         # Get 2d plane normal
# #
# #         # Get the 3d plane vecs centered on atom
# #
# #         # Get 3d plane normal
# #
# #
# #
# #     interpoland = np.array(
# #         [
# #             np.concatenate(
# #                 [
# #                     # tr[1].flatten(),
# #                     # tr[2].flatten(),
# #                     tr[3].flatten()
# #                 ]
# #             )
# #             for tr
# #             in transforms.values()
# #             ]
# #         )
# #     print(interpoland)
# #
# #     # transform_interpolator = LinearNDInterpolator(
# #     #     coord_array,
# #     #     interpoland
# #     # )
# #     transform_interpolator = RBFInterpolator(
# #         coord_array,
# #         interpoland,
# #         smoothing=2
# #     )
# #     # transform_interpolator_nearest = NearestNDInterpolator(
# #     #     coord_array,
# #     #     interpoland
# #     # )
# #
# #     # For each point
# #     values = []
# #     for sample in grid_samples:
# #         # for j, atom_id in enumerate(atom_ids):
# #         #     print(sample)
# #         # Get the anchor atoms
# #         dists, nbs = kd.query(sample, k=3)
# #         if dists[0] > 2:
# #             # print(dists[0])
# #             values.append(0)
# #             continue
# #         # print(nbs)
# #         # print(np.array(atom_ids)[nbs])
# #         # exit()
# #
# #         # Get the 2d Poss
# #         nbr_poss = {}
# #         for nbr in nbs:
# #             pos = coord_array[nbr]
# #             nbr_poss[atom_ids[nbr]] = pos
# #
# #         # Get the transform
# #         # tr = transforms[atom_ids[nbs[0]]]
# #         # tr_array = transform_interpolator(sample)[0]
# #         tr_array = transform_interpolator(sample.reshape(1,2))[0]
# #
# #         # print(f"Transform array")
# #         # print(tr_array)
# #         # if np.isnan(tr_array[0]):
# #         # tr_array = transform_interpolator_nearest(sample)[0]
# #             # print(tr_array)
# #
# #         # mat = tr_array[:9].reshape(3,3)
# #         # ref = tr_array[9:12]
# #         # mov = tr_array[12:16]
# #         # rot = tr_array[:3]
# #         # ref = tr_array[3:6]
# #         # mov = tr_array[6:9]
# #
# #         # print(mat)
# #         # print(ref)
# #         # print(mov)
# #
# #
# #         # Get the sample point
# #         # sample_point_2d = np.array([sample[0], sample[1], 0.0])
# #         # sample_point_2d_rel = sample_point_2d - ref
# #         # sample_point_3d_rel = np.matmul(
# #         #     spatial.transform.Rotation.from_euler('zyx', rot).as_matrix(),
# #         #     sample_point_2d_rel,
# #         # )
# #         # point_3d = sample_point_3d_rel + mov
# #         point_3d =
# #
# #         # Interpolate
# #         value = dmap.interpolate_value(
# #             # value=dmap.tricubic_interpolation(
# #             gemmi.Position(
# #                 point_3d[0],
# #                 point_3d[1],
# #                 point_3d[2],
# #             )
# #         )
# #         values.append(
# #             value
# #         )
# #
# #         if dists[0] < 1.0:
# #             rprint({
# #             "Pos": sample_point_2d,
# #             "Anchor Pos 2d": tr[2],
# #             "Anchor Pos 3d": tr[3],
# #             "Relative Pos": sample_point_2d_rel,
# #             "Relative Pos Distance": np.linalg.norm(sample_point_2d_rel),
# #             "Point 3D": point_3d,
# #             "Point 3d Relative": sample_point_3d_rel,
# #             "Point 3d Relative Dist": np.linalg.norm(sample_point_3d_rel),
# #             # "Reconstruction 2": np.dot(mat, components)
# #             "Value": value
# #         })
# #
# #
# #     plt.figure(figsize=(16, 9))
# #
# #     h = plt.scatter(
# #         grid_samples[:,0],
# #         grid_samples[:,1],
# #         c=values
# #     )
# #     plt.scatter(
# #         coord_array[:,0],
# #         coord_array[:,1],
# #         c='r'
# #     )
# #     plt.axis('scaled')
# #     plt.colorbar()
# #     print(f"Writing map!")
# #     plt.savefig('test.png')
#
# def plot_projection(structure_path,
#         cif_path,
#         map_path
# ):
#     # Load the structure
#     st = gemmi.read_structure(str(structure_path))
#
#     # Structure atom array
#     st_atom_pos_dict = get_structure_atom_array(st)
#
#     # Load the cif
#     cif = gemmi.cif.read(str(cif_path))
#     mol, atom_ids, atom_type_loop = cif_to_mol(cif)
#
#     # # Get the atom id array
#     # atom_id_array = get_atom_id_array(mol)
#
#     # Load the map
#     ccp4 = gemmi.read_ccp4_map(str(map_path), )
#     ccp4.setup(0.0)
#     dmap = ccp4.grid
#
#     # Generate the 2d projection
#     AllChem.Compute2DCoords(mol)
#     mask = []
#     for el in atom_type_loop:
#         if el == "H":
#             mask.append(False)
#         else:
#             mask.append(True)
#     coord_array = get_coord_array(mol)[np.array(mask), :]
#     print(coord_array)
#     print([np.min(coord_array, axis=0), np.max(coord_array, axis=0)])
#
#     # Get bounding box
#     bounds = get_bounds(coord_array)
#     print(f"Bounds: {bounds}")
#
#     # Generate the grid
#
#     xs = np.linspace(bounds[0][0], bounds[1][0], 100)
#     ys = np.linspace(bounds[0][1], bounds[1][1], 100)
#     grid_samples = np.array(
#         [x for x in itertools.product(
#             xs, ys
#         )]
#     )
#     print(grid_samples)
#     print(grid_samples.shape)
#     # Get the voronoi cells of the grid points relative to projection
#     # Get the atom coord kdtree
#     kd = KDTree(
#         coord_array
#     )
#
#     # For each atom get the three closest
#     transforms = {}
#     for j, pos in enumerate(coord_array):
#         dists, nbs = kd.query(pos, k=3)
#         poss_2d_unpadded = np.array(
#             [
#                 coord_array[nbr]
#                 for nbr
#                 in nbs
#                 ]
#             )
#
#         print(poss_2d_unpadded)
#
#         poss_2d = np.pad(
#             poss_2d_unpadded,
#             [(0, 0), (0, 1)]
#         )
#         poss_3d = np.array(
#             [
#             st_atom_pos_dict[atom_ids[nbr]]
#             for nbr
#             in nbs
#             ]
#         )
#         print(poss_2d)
#         print(poss_3d)
#         transform = get_transform(
#             poss_2d,
#             poss_3d
#         )
#         rprint(transform)
#         transforms[atom_ids[j]] = transform
#     for tr in transforms.values():
#         print(tr[1])
#
#     interpoland = np.array(
#         [
#             np.concatenate(
#                 [
#                     # tr[1].flatten(),
#                     # tr[2].flatten(),
#                     tr[3].flatten()
#                 ]
#             )
#             for tr
#             in transforms.values()
#             ]
#         )
#     print(interpoland)
#
#
#
#     transform_interpolator = LinearNDInterpolator(
#         coord_array,
#         interpoland
#     )
#     # transform_interpolator = RBFInterpolator(
#     #     coord_array,
#     #     interpoland,
#     #     # smoothing=2
#     # )
#     transform_interpolator_nearest = NearestNDInterpolator(
#         coord_array,
#         interpoland
#     )
#
#     # For each point
#     values = []
#     for sample in grid_samples:
#         # for j, atom_id in enumerate(atom_ids):
#         #     print(sample)
#         # Get the anchor atoms
#         dists, nbs = kd.query(sample, k=3)
#         if dists[0] > 2:
#             # print(dists[0])
#             values.append(0)
#             continue
#         # print(nbs)
#         # print(np.array(atom_ids)[nbs])
#         # exit()
#
#         # Get the 2d Poss
#         nbr_poss = {}
#         for nbr in nbs:
#             pos = st_atom_pos_dict[atom_ids[nbr]]
#             nbr_poss[atom_ids[nbr]] = pos
#
#         # Get the transform
#         # tr = transforms[atom_ids[nbs[0]]]
#         # tr_array = transform_interpolator(sample)[0]
#         tr_array = transform_interpolator(sample.reshape(1,2))[0]
#
#         # print(f"Transform array")
#         # print(tr_array)
#         if np.isnan(tr_array[0]):
#             tr_array = transform_interpolator_nearest(sample)[0]
#             # print(tr_array)
#
#         # mat = tr_array[:9].reshape(3,3)
#         # ref = tr_array[9:12]
#         # mov = tr_array[12:16]
#         # rot = tr_array[:3]
#         # ref = tr_array[3:6]
#         # mov = tr_array[6:9]
#
#         # print(mat)
#         # print(ref)
#         # print(mov)
#
#
#         # Get the sample point
#         # sample_point_2d = np.array([sample[0], sample[1], 0.0])
#         # sample_point_2d_rel = sample_point_2d - ref
#         # sample_point_3d_rel = np.matmul(
#         #     spatial.transform.Rotation.from_euler('zyx', rot).as_matrix(),
#         #     sample_point_2d_rel,
#         # )
#         # point_3d = sample_point_3d_rel + mov
#
#         point_3d = tr_array
#
#         # Interpolate
#         value = dmap.interpolate_value(
#             # value=dmap.tricubic_interpolation(
#             gemmi.Position(
#                 point_3d[0],
#                 point_3d[1],
#                 point_3d[2],
#             )
#         )
#         values.append(
#             value
#         )
#
#         if dists[0] < 0.5:
#             rprint({
#             "Pos": sample,
#                 "Anchor Pos 2d": coord_array[nbs[0]],
#             "Anchor Pos 3d": nbr_poss[atom_ids[nbs[0]]],
#             "Point 3D": point_3d,
#             "Value": value
#         })
#
#
#     plt.figure(figsize=(16, 9))
#
#     h = plt.scatter(
#         grid_samples[:,0],
#         grid_samples[:,1],
#         c=values
#     )
#     plt.scatter(
#         coord_array[:,0],
#         coord_array[:,1],
#         c='r'
#     )
#     plt.axis('scaled')
#     plt.colorbar()
#     print(f"Writing map!")
#     plt.savefig('test.png')

def plot_projection(structure_path,
        cif_path,
        map_path
):
    # Load the structure
    st = gemmi.read_structure(str(structure_path))

    # Structure atom array
    st_atom_pos_dict = get_structure_atom_array(st)

    # Load the cif
    cif = gemmi.cif.read(str(cif_path))
    mol, atom_ids, atom_type_loop = cif_to_mol(cif)

    # # Get the atom id array
    # atom_id_array = get_atom_id_array(mol)

    # Load the map
    ccp4 = gemmi.read_ccp4_map(str(map_path), )
    ccp4.setup(0.0)
    dmap = ccp4.grid

    # Generate the 2d projection
    AllChem.Compute2DCoords(mol)
    mask = []
    for el in atom_type_loop:
        if el == "H":
            mask.append(False)
        else:
            mask.append(True)
    coord_array = get_coord_array(mol)[np.array(mask), :]
    print(coord_array)
    print([np.min(coord_array, axis=0), np.max(coord_array, axis=0)])

    # Get bounding box
    bounds = get_bounds(coord_array)
    print(f"Bounds: {bounds}")

    # Generate the grid

    xs = np.linspace(bounds[0][0], bounds[1][0], 100)
    ys = np.linspace(bounds[0][1], bounds[1][1], 100)
    grid_samples = np.array(
        [x for x in itertools.product(
            xs, ys
        )]
    )
    print(grid_samples)
    print(grid_samples.shape)
    # Get the voronoi cells of the grid points relative to projection
    # Get the atom coord kdtree
    kd = KDTree(
        coord_array
    )

    # Get the 3d best fit plane thorugh pca
    pca_3d = decomposition.PCA(n_components=2)
    coord_array_3d = np.array(
        [
            x for x in st_atom_pos_dict.values()
        ]
    )
    pca_3d.fit(coord_array_3d)
    components_3d = pca_3d.components_
    mean_3d = np.mean(components_3d, axis=0)
    rprint(components_3d)
    plane_normal = np.cross(
        components_3d[0],
        components_3d[1]
    )
    rprint(plane_normal)


    # Get the
    pca_2d = decomposition.PCA(n_components=2)
    pca_2d.fit(coord_array)
    components_2d = pca_2d.components_
    mean_2d = np.mean(coord_array, axis=0)
    rprint(components_2d)

    # For each point
    values = []
    for sample in grid_samples:
        # Project point into ligand components
        sample_rel_2d  = pca_2d.transform((sample - mean_2d).reshape(1,-1))
        # rprint(f"Relative sample pos: {sample_rel_2d}")

        # Get the equivilent 3d pos
        sample_3d = pca_3d.inverse_transform(sample_rel_2d) + mean_3d
        # rprint(f"Sample pos 3d: {sample_3d}")
        # exit()

        # for j, atom_id in enumerate(atom_ids):
        #     print(sample)
        # Get the anchor atoms



        # Interpolate
        tmp_vals = []

        for x in np.linspace(-0.5,0.5,11):
            pos_to_sample = sample_3d[0] + (x*plane_normal)
            value = dmap.interpolate_value(
                # value=dmap.tricubic_interpolation(
                gemmi.Position(
                    pos_to_sample[0],
                    pos_to_sample[1],
                    pos_to_sample[2],
                )
            )
            tmp_vals.append(value)

        values.append(
            max(tmp_vals)
        )

        # if dists[0] < 0.5:
        #     rprint({
        #     "Pos": sample,
        #         "Anchor Pos 2d": coord_array[nbs[0]],
        #     "Anchor Pos 3d": nbr_poss[atom_ids[nbs[0]]],
        #     "Point 3D": point_3d,
        #     "Value": value
        # })


    a = plane_normal
    d = np.dot(mean_3d, plane_normal)
    projected_atom_poss = []
    for j, coord in enumerate(st_atom_pos_dict.values()):
        if atom_type_loop[j] == "H":
            continue
        plane_pos = coord - (
                (
            (np.dot(coord, a) - d) / np.dot(a,a)
        ) * a
        )
        projected_atom_poss.append(plane_pos)

    projected_atom_poss_array = np.array(projected_atom_poss)
    projected_atom_poss_array_rel = projected_atom_poss_array - mean_3d
    projected_atom_poss_array_rel_pca = pca_3d.transform(projected_atom_poss_array_rel)
    projected_atom_poss_array_2d = pca_2d.inverse_transform(projected_atom_poss_array_rel_pca) + mean_2d
    rprint(projected_atom_poss_array_2d)

    plt.figure(figsize=(16, 9))


    # h = plt.scatter(
    #     grid_samples[:,0],
    #     grid_samples[:,1],
    #     c=values
    # )
    plt.contourf(
        xs,
        ys,
        np.array(values).reshape(100, 100).T
    )
    # plt.imshow(np.array(values).reshape(100,100).T,
    #            extent=(
    #                bounds[0][0],
    #                bounds[1][0],
    #                bounds[0][1],
    #                bounds[1][1]))
    plt.scatter(
        projected_atom_poss_array_2d[:,0],
        projected_atom_poss_array_2d[:,1],
        c='r'
    )
    plt.axis('scaled')
    plt.colorbar()
    print(f"Writing map!")
    plt.savefig('test.png')