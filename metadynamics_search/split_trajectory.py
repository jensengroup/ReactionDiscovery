#!/home/mariaharris/.conda/envs/my-rdkit-env/bin/python

import os
import sys
import pickle
import xyz2mol_local
import reaction_box

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdmolops import GetFormalCharge
from reaction_box import reorder_atoms_to_map

def extract_smiles(xyz_file, charge):
    """
    uses xyz2mol to extract smiles with as much 3d structural information as
    possible
    """
    atoms, _, xyz_coordinates = xyz2mol_local.read_xyz_file(xyz_file)
    input_mol = xyz2mol_local.xyz2mol(atoms, xyz_coordinates, charge=charge,
                                      use_graph=True,
                                      allow_charged_fragments=True,
                                      use_huckel=True, use_atom_maps=True,
                                      embed_chiral=True)
    input_mol = reorder_atoms_to_map(input_mol)
    structure_mol = reaction_box.choose_resonance_structure(input_mol)
    structure_mol = reaction_box.chiral_tags(structure_mol)
    rdmolops.AssignStereochemistry(structure_mol)
    structure_smiles = Chem.MolToSmiles(structure_mol)

    return structure_smiles, GetFormalCharge(structure_mol)

def split_trajectory(trajectory_file, n_check, charge):
    """
    Check trajectory for reactions every n_check points
    """

    os.mkdir('analyze_trajectory')
    smiles_list = []
    n_steps_list = []
    count = 0
    saved_struc = 1
    n_steps = 0
    with open(trajectory_file, 'r') as _file:
        line = _file.readline()
        n_lines = int(line.split()[0])+2
        while line:
            if count % (n_lines*n_check) == 0:
                file_name = 'analyze_trajectory/'+str(saved_struc)+'.xyz'
                with open(file_name, 'w') as xyz_file:
                    for _ in range(n_lines):
                        xyz_file.write(line)
                        line = _file.readline()
                        count += 1
                try:
                    smiles, formal_charge = extract_smiles(file_name, charge)
                    if not smiles_list and formal_charge == charge:
                        smiles_list.append(smiles)
                        print('smiles =', smiles)
                        saved_struc += 1
                        n_steps_list.append(n_steps)
                    elif smiles != smiles_list[-1] and formal_charge == charge:
                        smiles_list.append(smiles)
                        print('smiles = ', smiles)
                        saved_struc += 1
                        n_steps_list.append(n_steps)
                except (ValueError, RuntimeError):
                    pass
                n_steps += n_check
            else:
                line = _file.readline()
                count += 1
    print(n_steps_list)
    return n_steps_list

if __name__ == '__main__':
    trajectory_file = sys.argv[1]
    n_check = int(sys.argv[2])
    charge = int(sys.argv[3])
    outfile = sys.argv[4]
    n_steps_list = split_trajectory(trajectory_file, n_check, charge)
    with open(outfile, 'wb') as _file:
        pickle.dump(n_steps_list, _file)

