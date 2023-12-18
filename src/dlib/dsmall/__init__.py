from pathlib import Path

import numpy as np
import gemmi
import networkx
from rdkit import Chem

from .. import constants
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
        nearest_image = st1.cell.find_nearest_image(atom_1.pos, atom_2.pos)
        distance = nearest_image.dist()

        # distance = atom_1.pos.dist(atom_2.pos)
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
            nearest_image = st1.cell.find_nearest_image(atom1.pos, atom2.pos)
            distance = nearest_image.dist()
            # atom_distances.append(atom1.pos.dist(atom2.pos))
            atom_distances.append(distance)
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
        if atom_1.element.name == "H":
            continue

        match = None

        for atom_2 in res_2:
            if atom_1.name == atom_2.name:
                match = atom_2.name
        if not match:
            return None
        atom_matches.append(
            (atom_1.name, match)
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


def get_fragment_mol_from_dataset_cif_path(dataset_cif_path: Path):
    # Open the cif document with gemmi
    cif = gemmi.cif.read(str(dataset_cif_path))

    # Create a blank rdkit mol
    mol = Chem.Mol()
    editable_mol = Chem.EditableMol(mol)

    key = "comp_LIG"
    try:
        cif['comp_LIG']
    except:
        key = "data_comp_XXX"

    # Find the relevant atoms loop
    atom_id_loop = list(cif[key].find_loop('_chem_comp_atom.atom_id'))
    atom_type_loop = list(cif[key].find_loop('_chem_comp_atom.type_symbol'))
    atom_charge_loop = list(cif[key].find_loop('_chem_comp_atom.charge'))
    if not atom_charge_loop:
        atom_charge_loop = list(cif[key].find_loop('_chem_comp_atom.partial_charge'))
        if not atom_charge_loop:
            atom_charge_loop = [0]*len(atom_id_loop)

    aromatic_atom_loop = list(cif[key].find_loop('_chem_comp_atom.aromatic'))
    if not aromatic_atom_loop:
        aromatic_atom_loop = [None]*len(atom_id_loop)

    # Get the mapping
    id_to_idx = {}
    for j, atom_id in enumerate(atom_id_loop):
        id_to_idx[atom_id] = j

    # Iteratively add the relveant atoms
    for atom_id, atom_type, atom_charge in zip(atom_id_loop, atom_type_loop, atom_charge_loop):
        if len(atom_type) > 1:
            atom_type = atom_type[0] + atom_type[1].lower()
        atom = Chem.Atom(atom_type)
        atom.SetFormalCharge(round(float(atom_charge)))
        editable_mol.AddAtom(atom)

    # Find the bonds loop
    bond_1_id_loop = list(cif[key].find_loop('_chem_comp_bond.atom_id_1'))
    bond_2_id_loop = list(cif[key].find_loop('_chem_comp_bond.atom_id_2'))
    bond_type_loop = list(cif[key].find_loop('_chem_comp_bond.type'))
    aromatic_bond_loop = list(cif[key].find_loop('_chem_comp_bond.aromatic'))
    if not aromatic_bond_loop:
        aromatic_bond_loop = [None]*len(bond_1_id_loop)

    try:
        # Iteratively add the relevant bonds
        for bond_atom_1, bond_atom_2, bond_type, aromatic in zip(bond_1_id_loop, bond_2_id_loop, bond_type_loop, aromatic_bond_loop):
            bond_type = constants.BondTypeCifToRdkit[bond_type]
            if aromatic:
                if aromatic == "y":
                    bond_type = constants.BondTypeCifToRdkit['aromatic']

            editable_mol.AddBond(
                id_to_idx[bond_atom_1],
                id_to_idx[bond_atom_2],
                order=bond_type
            )
    except Exception as e:
        print(e)
        print(atom_id_loop)
        print(id_to_idx)
        print(bond_1_id_loop)
        print(bond_2_id_loop)
        raise Exception

    edited_mol = editable_mol.GetMol()
    # for atom in edited_mol.GetAtoms():
    #     print(atom.GetSymbol())
    #     for bond in atom.GetBonds():
    #         print(f"\t\t{bond.GetBondType()}")
    # for bond in edited_mol.GetBonds():
    #     ba1 = bond.GetBeginAtomIdx()
    #     ba2 = bond.GetEndAtomIdx()
    #     print(f"{bond.GetBondType()} : {edited_mol.GetAtomWithIdx(ba1).GetSymbol()} : {edited_mol.GetAtomWithIdx(ba2).GetSymbol()}")  #*}")
    # print(Chem.MolToMolBlock(edited_mol))


    # HANDLE SULFONATES
    # forward_mol = Chem.ReplaceSubstructs(
    #     edited_mol,
    #     Chem.MolFromSmiles('S(O)(O)(O)'),
    #     Chem.MolFromSmiles('S(=O)(=O)(O)'),
    #     replaceAll=True,)[0]
    patt = Chem.MolFromSmarts('S(-O)(-O)(-O)')
    matches = edited_mol.GetSubstructMatches(patt)

    sulfonates = {}
    for match in matches:
        sfn = 1
        sulfonates[sfn] = {}
        on = 1
        for atom_idx in match:
            atom = edited_mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "S":
                sulfonates[sfn]["S"] = atom_idx
            else:
                atom_charge = atom.GetFormalCharge()

                if atom_charge == -1:
                    continue
                else:
                    if on == 1:
                        sulfonates[sfn]["O1"] = atom_idx
                        on += 1
                    elif on == 2:
                        sulfonates[sfn]["O2"] = atom_idx
                        on += 1
                # elif on == 3:
                #     sulfonates[sfn]["O3"] = atom_idx
    print(f"Matches to sulfonates: {matches}")

    # atoms_to_charge = [
    #     sulfonate["O3"] for sulfonate in sulfonates.values()
    # ]
    # print(f"Atom idxs to charge: {atoms_to_charge}")
    bonds_to_double =[
        (sulfonate["S"], sulfonate["O1"]) for sulfonate in sulfonates.values()
    ] + [
        (sulfonate["S"], sulfonate["O2"]) for sulfonate in sulfonates.values()
    ]
    print(f"Bonds to double: {bonds_to_double}")

    # Replace the bonds and update O3's charge
    new_editable_mol = Chem.EditableMol(Chem.Mol())
    for atom in edited_mol.GetAtoms():
        atom_idx = atom.GetIdx()
        new_atom = Chem.Atom(atom.GetSymbol())
        charge = atom.GetFormalCharge()
        # if atom_idx in atoms_to_charge:
        #     charge = -1
        new_atom.SetFormalCharge(charge)
        new_editable_mol.AddAtom(new_atom)

    for bond in edited_mol.GetBonds():
        bond_atom_1 = bond.GetBeginAtomIdx()
        bond_atom_2 = bond.GetEndAtomIdx()
        double_bond = False
        for bond_idxs in bonds_to_double:
            if (bond_atom_1 in bond_idxs) & (bond_atom_2 in bond_idxs):
                double_bond = True
        if double_bond:
            new_editable_mol.AddBond(
                bond_atom_1,
                bond_atom_2,
                order=constants.BondTypeCifToRdkit['double']
            )
        else:
            new_editable_mol.AddBond(
                bond_atom_1,
                bond_atom_2,
                order=bond.GetBondType()
            )
    new_mol = new_editable_mol.GetMol()
    # print(Chem.MolToMolBlock(new_mol))

    new_mol.UpdatePropertyCache()
    # Chem.SanitizeMol(new_mol)
    return new_mol