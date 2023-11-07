import numpy as np
import gemmi
import networkx

from ..dcommon import AtomID

AtomMatch = list[tuple[AtomID, AtomID]]

def get_structure_atom(st, atom_id):
    model = st[0]
    chain = model[atom_id[0]]
    res_group = chain[atom_id[1]]
    res = res_group[0]
    atom_group = res[atom_id[2]]
    atom = atom_group[0]
    return atom

def get_structure_res(st, residue_id):
    model = st[0]
    chain = model[residue_id[0]]
    res_group = chain[residue_id[1]]
    res = res_group[0]
    return res

def get_rmsd_from_match(
        st1,
        st2,
        match
):
    distances = []
    for atom_1_id, atom_2_id in match:
        atom_1 = get_structure_atom(st1, atom_1_id)
        atom_2 = get_structure_atom(st2, atom_2_id)
        distance = atom_1.pos.dist(atom_2.pos)
        distances.append(distance)

    return np.sqrt(np.mean(np.square(distances)))

def get_rmsd_from_closest_atom(
            st1,
            st2,
            st1_lig_id,
            st2_lig_id
        ):

    lig_1_res = get_structure_res(st1, st1_lig_id)
    lig_2_res = get_structure_res(st2, st2_lig_id)

    distances = []
    for atom1 in lig_1_res:
        atom_distances = []
        for atom2 in lig_2_res:
            atom_distances.append(atom1.pos.dist(atom2.pos))
        distances.append(min(atom_distances))

    return np.sqrt(np.mean(np.square(distances)))


def get_ligands(st):
    ligands = {}
    for model in st:
        for chain in model:
            for res in chain:
                if res.name == "LIG":
                    ligands[(chain.name, str(res.seqid.num))] = res

    return ligands

def cif_to_graph(cif):
    key = "comp_LIG"
    try:
        cif['comp_LIG']
    except:
        key = "data_comp_XXX"

    atom_id_loop = list(cif[key].find_loop('_chem_comp_atom.atom_id'))
    atom_type_loop = list(cif[key].find_loop('_chem_comp_atom.type_symbol'))
    bond_1_id_loop = list(cif[key].find_loop('_chem_comp_bond.atom_id_1'))
    bond_2_id_loop = list(cif[key].find_loop('_chem_comp_bond.atom_id_2'))

    graph = networkx.Graph()

    for atom_name, atom_type in zip(atom_id_loop, atom_type_loop):
        graph.add_node(atom_name, Z=atom_type)

    for bond_1_atom_id, bond_2_atom_id in zip(bond_1_id_loop, bond_2_id_loop):
        graph.add_edge(bond_1_atom_id, bond_2_atom_id)  # ignoring bond type

    return graph

def match_res_to_cif(res, cif):
    ...

def get_match(res_1, res_2):
    atom_matches = []
    for atom_1 in res_1:
        match = None

        for atom_2 in res_2:
            if atom_1.name == atom_2.name:
                match = atom_2.name
        if not match:
            return None
        atom_matches.append(
            (atom_1.name, atom_2.name)
        )

    return atom_matches


def match_structure_ligands(
        structure_1,
        structure_2,
        # cif
):
    # Get structure 1 ligands
    st1_ligands = get_ligands(structure_1)

    # Get structure 2 ligands
    st2_ligands = get_ligands(structure_2)

    # # Get the cif graph
    # graph = cif_to_graph(cif)

    # Try all matches
    matches = {}
    for st1_ligand_id, st1_ligand in st1_ligands.items():

        # st1_cif_match: AtomMatch = match_res_to_cif(st1_ligand, cif)

        for st2_ligand_id, st2_ligand in st2_ligands.items():

            # st2_cif_match: AtomMatch = match_res_to_cif(st2_ligand, cif)

            matched_atoms = get_match(
                st1_ligand,
                st2_ligand
            )
            if not matched_atoms:
                continue
            # for atom_pair in matched_atoms:
            matches[(st1_ligand_id, st2_ligand_id)] = [
                (
                    (st1_ligand_id[0], st1_ligand_id[1], atom_pair[0],),
                    (st2_ligand_id[0], st2_ligand_id[1], atom_pair[1],),
                )
                for atom_pair
                in matched_atoms
            ]

    return matches