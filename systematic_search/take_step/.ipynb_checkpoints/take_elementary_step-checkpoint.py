from itertools import product, combinations
import numpy as np
import copy
from tqdm.notebook import tqdm

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
from rdkit.Chem.rdchem import ResonanceFlags
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

from .compound import Compound

def valid_Is(AC, atomic_num, Cs):
    """ Tests is atoms have larger valence than allowed max valence """

    max_valence = {}

    max_valence[1] = 1
    max_valence[6] = 4
    max_valence[7] = 4
    max_valence[8] = 3  
    max_valence[9] = 1
    max_valence[14] = 4
    max_valence[15] = 5
    max_valence[16] = 6
    max_valence[17] = 1
    max_valence[35] = 1
    max_valence[53] = 1
    
    # right now this is not pretty
    Is = []
    for C in Cs:
        I = AC + C
        ok = True
        for atom, valence in zip(atomic_num, I.sum(axis=1, dtype=np.int)):
            if valence > max_valence[atom]:
                ok = False
        if ok:
            Is.append(I)
                    
    return Is


def get_valid_Is(AC, atomic_num, max_bonds=2):
    """ Get valid I matrices not optimized on memory or speed but works
    TODO: use cython to speed up, and sparse matrices to reduce mem.
    """
    num_atoms = len(atomic_num)
    Cs = []
    make1, break1 = [], []
    for i in range(num_atoms):
        for j in range(num_atoms):
            C = np.zeros((num_atoms,num_atoms), dtype=np.int32) # make sparse matrix here
            if j > i:
                if AC[i,j] == 0:
                    C[i,j] = C[j,i] = 1
                    make1.append(C)
                else:
                     C[i,j] = C[j,i] = -1
                     break1.append(C)
    
    make1_break1 = [x[0] + x[1] for x in product(make1, break1)]
    Cs += make1 + break1 + make1_break1

    if max_bonds >=2:
        make2 = [x[0] + x[1] for x in combinations(make1, 2)]
        break2 = [x[0] + x[1] for x in combinations(break1, 2)]
        make2_break1 = [x[0] + x[1] for x in product(make2, break1)]
        make1_break2 = [x[0] + x[1] for x in product(make1, break2)]
        make2_break2 = [x[0] + x[1] for x in product(make2, break2)]
        Cs += make2 + break2 + make2_break1 + make1_break2 + make2_break2
    
    if max_bonds >= 3:
        make3 = [x[0] + x[1] + x[2] for x in combinations(make1, 3)]
        break3 = [x[0] + x[1] + x[2] for x in combinations(break1, 3)]
        make3_break1 = [x[0] + x[1] for x in product(make3, break1)]
        make3_break2 = [x[0] + x[1] for x in product(make3, break2)]
        make2_break3 = [x[0] + x[1] for x in product(make2, break3)]
        make1_break3 = [x[0] + x[1] for x in product(make1, break3)]        
        Cs += make3 + break3 + make3_break1 + make3_break2 + make2_break3 + make1_break3
    
    Is = valid_Is(AC, atomic_num, Cs)
    print(f">> Valid Is {len(Is)}")
    
    return Is


def reassign_atom_idx(mol):
    """ Assigns RDKit mol atomid to atom mapped id """ 
    renumber = [(atom.GetIdx(), atom.GetAtomMapNum()) for atom in mol.GetAtoms()]
    new_idx = [idx[0] for idx in sorted(renumber, key = lambda x: x[1])]
    
    return Chem.RenumberAtoms(mol, new_idx)


def set_chirality(product, reactant):
    """ Produce all combinations of isomers (R/S, and cis/trans). """ 

    # TODO move these somewhere it makes more sense.
    product = reassign_atom_idx(product)
    reactant = reassign_atom_idx(reactant)

    Chem.SanitizeMol(product)
    Chem.SanitizeMol(reactant)

    # find unchanged atoms (they need the same chirality) 
    # TODO: Something is wrong here.
    chiral_atoms_product = Chem.FindMolChiralCenters(product, includeUnassigned=True)

    patt = '[C^3;H2,H1]'
    sp3_hyp_C = Chem.MolFromSmarts(patt)

    mcs_product = product.GetSubstructMatch(sp3_hyp_C)
    mcs_reactant = reactant.GetSubstructMatch(sp3_hyp_C)

    unchanged_atoms = []
    for atom, chiral_tag in chiral_atoms_product:
        product_neighbors = [a.GetIdx() for a in product.GetAtomWithIdx(atom).GetNeighbors()]
        reactant_neighbors = [a.GetIdx() for a in reactant.GetAtomWithIdx(atom).GetNeighbors()]
        
        if sorted(product_neighbors) == sorted(reactant_neighbors):
            unchanged_atoms.append(atom)

    # make combinations of isomers.
    opts = StereoEnumerationOptions(onlyUnassigned=False, unique=False) # should i use unique?
    rdmolops.AssignStereochemistry(product, cleanIt=True, flagPossibleStereoCenters=True, force=True)

    product_isomers = []
    product_isomers_mols = []
    for product_isomer in EnumerateStereoisomers(product, options=opts):
        rdmolops.AssignStereochemistry(product_isomer, force=True)
        for atom in unchanged_atoms:
            reactant_global_tag = reactant.GetAtomWithIdx(atom).GetProp('_CIPCode') # R/S ?

            # TODO make sure that the _CIPRank is the same for atom in reactant and product.
            
            product_isomer_global_tag = product_isomer.GetAtomWithIdx(atom).GetProp('_CIPCode') # R/S ? 
            if reactant_global_tag != product_isomer_global_tag:
                product_isomer.GetAtomWithIdx(atom).InvertChirality()

        if Chem.MolToSmiles(product_isomer) not in product_isomers:
            product_isomers.append(Chem.MolToSmiles(product_isomer))
            product_isomers_mols.append(product_isomer)

    return product_isomers_mols


def take_elementary_step(reac_compound, max_bonds=2, use_atom_maps=True, all_isomers=True, 
            all_resonance=False):
    """ Takes elementary step. This function returns all possible 
    combinations, that have a valled valence """

    # 1) find the ac matrix of the reactant
    reactant_mol = reac_compound.get_rdkit_mol()
    ac_matrix = Chem.GetAdjacencyMatrix(reactant_mol, useBO=False)
    atomic_num_list = Compound.symbols2numbers(reac_compound.atomic_symbols)
  
    # 2) build all possible conversion matrices from ac_matrix, allowing max_bonds bonds 
    # to break and form. 
    # TODO: rewrite to cython/c++ and use sparse matrices. And eleminate non-valid products 
    # when they are made - saves memory. Or åerhaps just as a generator.
    I_elementary = get_valid_Is(ac_matrix, atomic_num_list, max_bonds=max_bonds) 
    
    all_compounds = []
    all_smiles = []
    count = 0
    for I in tqdm(I_elementary, total=len(I_elementary)):
        new_mol = Compound.from_ac_matrix(I, reac_compound.atomic_symbols, 
                                            addHs=False, charge=reac_compound.charge,
                                            multiplicity=reac_compound.multiplicity,
                                            charged_fragments=True, 
                                            use_atom_maps=use_atom_maps)

        # This is a hack for Grambow mol
        try:
            Chem.SanitizeMol(new_mol.rdkit_mol)
        except:
            continue
            
        # Check that the total charge is 0. Hack for Grambow.
        if int(Chem.GetFormalCharge(new_mol.rdkit_mol)) != 0:
            continue
            
        # make all resonance forms:
        all_resonance_structures = [res for res in rdchem.ResonanceMolSupplier(new_mol.rdkit_mol, 
                        ResonanceFlags.UNCONSTRAINED_ANIONS)] # ,  ResonanceFlags.ALLOW_CHARGE_SEPARATION
        if all_resonance: # Hack to compare with Grambow.
            resonance_structures = all_resonance_structures
        else:
            if len(all_resonance_structures) == 1:
                resonance_structures = all_resonance_structures
            else:
                most_rigid_res, min_rot_bonds = None, 999
                for res in all_resonance_structures:
                    Chem.SanitizeMol(res)
                    num_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(res)
                    if num_rot_bonds < min_rot_bonds:
                        most_rigid_res = copy.deepcopy(res)
                        min_rot_bonds = num_rot_bonds
                resonance_structures = [most_rigid_res]

        # Assign Stero Chemistry - This will return a list of molecules.
        isomers = []
        if all_isomers:
            for resonance_structure in resonance_structures:
                isomers += set_chirality(resonance_structure, reactant_mol)
        else:
            isomers += [resonance_structures]

        for isomer in isomers: # consider not writing smiles two times.
            if Chem.MolToSmiles(isomer) not in all_smiles:
                all_smiles.append(Chem.MolToSmiles(isomer))
                isomer = Chem.MolFromMolBlock(Chem.MolToMolBlock(isomer),sanitize=False)
                Chem.SanitizeMol(isomer) # or update property chache
                all_compounds.append(isomer)

    return all_smiles, all_compounds 