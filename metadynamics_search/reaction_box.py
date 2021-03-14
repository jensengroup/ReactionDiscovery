#!/home/mariaharris/.conda/envs/my-rdkit-env/bin/python

"""
Find reactant products for the input SMILES (atom-mapped) string using xTB meta-dynamics
"""

import os
import shutil
import time
import tarfile
import itertools
import pickle
import xyz2mol_local
import numpy as np
import pandas as pd
from rdkit import RDLogger
from rdkit.Chem.rdmolops import GetFormalCharge
from rdkit.Chem import rdmolops

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

from rdkit.Chem.rdmolops import FastFindRings

def reorder_atoms_to_map(mol):

    """
    Reorders the atoms in a mol objective to match that of the mapping
    """

    atom_map_order = np.zeros(mol.GetNumAtoms()).astype(int)
    for atom in mol.GetAtoms():
        map_number = atom.GetAtomMapNum()-1
        atom_map_order[map_number] = atom.GetIdx()
    mol = Chem.RenumberAtoms(mol, atom_map_order.tolist())
    return mol

def choose_resonance_structure(mol):
    """
    This function creates all resonance structures of the mol object, counts
    the number of rotatable bonds for each structure and chooses the one with
    fewest rotatable bonds (most 'locked' structure)
    """
    resonance_mols = rdchem.ResonanceMolSupplier(mol,
                                                 rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION)
    new_mol = None
    for res_mol in resonance_mols:
        Chem.SanitizeMol(res_mol)
        n_rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(res_mol)
        if new_mol is None:
            smallest_rot_bonds = n_rot_bonds
            new_mol = res_mol
        if n_rot_bonds < smallest_rot_bonds:
            smallest_rot_bonds = n_rot_bonds
            new_mol = res_mol

    Chem.DetectBondStereochemistry(new_mol, -1)
    rdmolops.AssignStereochemistry(new_mol, flagPossibleStereoCenters=True,
                                   force=True)
    Chem.AssignAtomChiralTagsFromStructure(new_mol, -1)
    return new_mol


def chiral_tags(mol):
    """
    Tag methylene and methyl groups with a chiral tag priority defined
    from the atom index of the hydrogens
    """
    li_list = []
    smarts_ch2 = '[!#1][*]([#1])([#1])([!#1])'
    atom_sets = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts_ch2))
    for atoms in atom_sets:
        atoms = sorted(atoms[2:4])
        prioritized_H = atoms[-1]
        li_list.append(prioritized_H)
        mol.GetAtoms()[prioritized_H].SetAtomicNum(9)
    smarts_ch3 = '[!#1][*]([#1])([#1])([#1])'
    atom_sets = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts_ch3))
    for atoms in atom_sets:
        atoms = sorted(atoms[2:])
        H1 = atoms[-1]
        H2 = atoms[-2]
        li_list.append(H1)
        li_list.append(H2)
        mol.GetAtoms()[H1].SetAtomicNum(9)
        mol.GetAtoms()[H2].SetAtomicNum(9)

    Chem.AssignAtomChiralTagsFromStructure(mol, -1)
    rdmolops.AssignStereochemistry(mol)
    for atom_idx in li_list:
        mol.GetAtoms()[atom_idx].SetAtomicNum(1)

    return mol


def write_xyz_file(mol, file_name):

    """
    Embeds a mol object to get 3D coordinates which are written to an .xyz file
    """

    n_atoms = mol.GetNumAtoms()
    charge = Chem.GetFormalCharge(mol)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]

    Chem.SanitizeMol(mol)
    rdmolops.AssignStereochemistry(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    with open(file_name, 'w') as _file:
        _file.write(str(n_atoms)+'\n\n')
        for atom, symbol in enumerate(symbols):
            coord = mol.GetConformers()[0].GetAtomPosition(atom)
            line = " ".join((symbol, str(coord.x), str(coord.y), str(coord.z), "\n"))
            _file.write(line)
        if charge != 0:
            _file.write("$set\n")
            _file.write("chrg "+str(charge)+"\n")
            _file.write("$end")


def submit_md_run(file_name, box_scale, reactant_idx):
    """
    submits a xTB md run and extracts last structure
    """
    jobid = os.popen('../submit_scripts/submit_xtb_md ' + file_name\
            + ' ' + str(box_scale)).read()
    jobid = int(jobid.split()[-1])
    return {jobid}

def check_finished_jobs(job_ids, job_dictionary):
    """
    for the finished jobs: check whether 
        - the job excited normally: remove jobid from checked tuple
        - the job excited abnormally: submit the job again
    """

    successful_jobs = set()
    crashed_run_numbers = []
    for jobid in job_ids:
        job_crashed = False
        run_number = job_dictionary[str(jobid)]
        os.chdir('run'+str(run_number))
        with open(str(run_number)+'.out', 'r') as ofile:
            line = ofile.readline()
            while line:
                if 'emergency exit' in line:
                    job_crashed = True
                    crashed_run_numbers.append(run_number)
                line = ofile.readline()
        os.chdir('../')
        if not job_crashed:
            successful_jobs = successful_jobs|{jobid}

    return successful_jobs, crashed_run_numbers


def wait_for_jobs_to_finish(job_ids):
    """
    This function checks with slurm if a specific set of jobids is finished with a
    frequency of 1 minute.
    Stops when the jobs are done.
    """
    while True:
        job_info1 = os.popen("squeue -p mko").readlines()[1:]
        job_info2 = os.popen("squeue -u mariaharris").readlines()[1:]
        current_jobs1 = {int(job.split()[0]) for job in job_info1}
        current_jobs2 = {int(job.split()[0]) for job in job_info2}
        current_jobs = current_jobs1|current_jobs2
        if current_jobs.isdisjoint(job_ids):
            break
        else:
            time.sleep(30)


def wait_for_jobs_to_finish_meta(job_ids, job_dictionary, xyz_file, kpush, alp,
                                 scale_factor):
    """
    This function checks with slurm if a specific set of jobids is finished with a
    frequency of 30 s.
    Stops when the jobs are done.
    """
    while True:
        job_info1 = os.popen("squeue -p mko").readlines()[1:]
        job_info2 = os.popen("squeue -u mariaharris").readlines()[1:]
        current_jobs1 = {int(job.split()[0]) for job in job_info1}
        current_jobs2 = {int(job.split()[0]) for job in job_info2}
        current_jobs = current_jobs1|current_jobs2

        if current_jobs.isdisjoint(job_ids):
            break

        intersection_jobs = current_jobs.intersection(job_ids)
        finished_jobs = job_ids.symmetric_difference(intersection_jobs)
        if finished_jobs:
            successful_jobs, crashed_run_numbers = check_finished_jobs(finished_jobs,
                                                                       job_dictionary)
            job_ids = job_ids.symmetric_difference(finished_jobs)
            for i in crashed_run_numbers:
                shutil.rmtree('run'+str(i))
            new_job_ids, new_job_dictionary = submit_metadynamics_runs(xyz_file, kpush,
                                                                       alp, scale_factor,
                                                                       crashed_run_numbers)
            job_ids = job_ids|new_job_ids
            job_dictionary.update(new_job_dictionary)

        time.sleep(60)

def extract_last_structure(trj_file, last_structure_name):
    """
    Extracts the last structure in a trajectory file
    """

    with open(trj_file, 'r') as _file:
        line = _file.readline()
        n_lines = int(line.split()[0])+2
    count = 0
    input_file = open(trj_file, 'r')
    dest = None
    for line in input_file:
        if count % n_lines == 0:
            if dest:
                dest.close()
            dest = open(last_structure_name, "w")
        count += 1
        dest.write(line)

def calculate_md_relaxed_structure(xyz_file, scale_factor, reactant_index):
    """
    This function submits an md with a box with size scaled by scale_factor and
    extracts last structure of the trajectory file
    """
    jobid = submit_md_run(xyz_file, scale_factor, reactant_index)
    wait_for_jobs_to_finish(jobid)
    trj_file = str(scale_factor)+'_md.trj'
    out_file = str(scale_factor)+'_md.xyz'
    extract_last_structure(trj_file, out_file)
    shutil.copy(out_file, '../')
    return out_file


def submit_metadynamics_runs(xyz_file, k_push, alp, box_scale, run_numbers):
    """
    submits a metadynamics calculation with the indicated parameters
    """
    job_dictionary={}
    jobids = set()
    for i in run_numbers:
        os.mkdir('run'+str(i))
        shutil.copy('metadyn.inp', 'run'+str(i))
        shutil.copy(xyz_file, 'run'+str(i))
        os.chdir('run'+str(i))
        jobid = os.popen('../../../submit_scripts/submit_xtb_meta ' + xyz_file + ' ' + str(k_push) \
        + ' ' + str(alp) + ' ' +  str(box_scale) + ' ' + str(i)).read()
        job_dictionary[jobid.split()[-1]] = str(i)
        jobid = {int(jobid.split()[-1])}
        jobids = jobids|jobid
        os.chdir('../')
    return jobids, job_dictionary

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
    structure_mol = choose_resonance_structure(input_mol)
    structure_mol = chiral_tags(structure_mol)
    rdmolops.AssignStereochemistry(structure_mol)
    structure_smiles = Chem.MolToSmiles(structure_mol)

    return structure_smiles, GetFormalCharge(structure_mol)


def canonicalize_smiles(structure_smiles):
    """
    remove all structural info an atom mapping information
    """
    mol = Chem.MolFromSmiles(structure_smiles, sanitize=False)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    canonical_smiles = Chem.MolToSmiles(mol)

    return canonical_smiles

def check_opt_xyz(charge, n_steps_list):

    n_steps_list_opt = []
    reactions = []
    canonical_reactants = []
    canonical_products = []
    smiles_list = []

    opt_files = [f for f in os.listdir(os.curdir) if f.endswith("optxyz")]
    opt_files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    for i, optfile in enumerate(opt_files[:-1]):
        try:
            smiles, formal_charge = extract_smiles(optfile, charge)
            print(smiles)
            if not smiles_list and formal_charge == charge:
                smiles_list.append(smiles)
                shutil.copy(optfile, 'optimized_structures')

            elif smiles != smiles_list[-1] and formal_charge == charge:
                reactant = smiles_list[-1]
                reaction = reactant+'>>'+smiles
                reactions.append(reaction)

                canonical_reactant = canonicalize_smiles(smiles_list[-1])
                canonical_reactants.append(canonical_reactant)
                canonical_product = canonicalize_smiles(smiles)
                canonical_products.append(canonical_product)

                smiles_list.append(smiles)
                n_steps_list_opt.append(n_steps_list[i])
                shutil.copy(optfile, 'optimized_structures')
        except ValueError:
            pass
    #print(n_steps_list_opt)


    return reactions, canonical_reactants, canonical_products, n_steps_list_opt



def check_metadyn_runs(n_runs, n_check, charge, combined_df_name):

    jobids = set()
    for i in range(n_runs):
        os.chdir('run'+str(i))
        jobid = os.popen('../../../submit_scripts/submit_python '+\
            '/home/mariaharris/github/ReactionDiscovery/metadynamics_search/split_trajectory.py xtb.trj '\
                         +str(n_check) + ' ' + str(charge) + ' list_file.p').read()
        jobid = {int(jobid.split()[-1])}
        jobids = jobids|jobid
        os.chdir('../')

    wait_for_jobs_to_finish(jobids)

    jobids = set()
    for i in range(n_runs):
        os.chdir('run'+str(i)+'/analyze_trajectory')
        jobid = os.popen('../../../../submit_scripts/submit_batches_xtb_opt').read()
        jobid = {int(jobid.split()[-1])}
        jobids = jobids|jobid
        os.chdir('../../')

    wait_for_jobs_to_finish(jobids)

    combined_df = pd.DataFrame()

    for i in range(n_runs):
        with open('run'+str(i)+'/list_file.p', 'rb') as list_file:
            n_steps_list = pickle.load(list_file)
        os.chdir('run'+str(i)+'/analyze_trajectory')
        tfs = [f for f in os.listdir(os.curdir) if f.endswith("log_xtb.tar.gz")]
        for filename in tfs:
            tf = tarfile.open(filename, mode='r')
            tf.extractall()
            tf.close()
        os.mkdir('optimized_structures')
        reactions, canonical_reactants, canonical_products, n_steps_list_opt = \
                check_opt_xyz(charge, n_steps_list)
        df = pd.DataFrame()
        df['Reactions'] = reactions
        df['Reactants_canonical'] = canonical_reactants
        df['Products_canonical'] = canonical_products
        df['N_steps'] = n_steps_list_opt
        df.to_pickle('dataframe.pkl')
        combined_df = pd.concat([combined_df, df])
        os.chdir('../../')

    combined_df_final = combined_df.groupby(['Reactions', 'Reactants_canonical',
                                             'Products_canonical']).size().reset_index(name='counts')

    print(combined_df_final)
    print(combined_df.groupby(['Reactions', 'Reactants_canonical',
                               'Products_canonical'])['N_steps'].mean())

    combined_df_final['average steps'] = combined_df.groupby(['Reactions', 'Reactants_canonical',
                                                              'Products_canonical'])['N_steps'].mean().values
    combined_df_final['median steps'] = combined_df.groupby(['Reactions', 'Reactants_canonical',
                                                             'Products_canonical'])['N_steps'].median().values
    print(combined_df_final)
    combined_df_final.to_pickle(combined_df_name)




def test_parameters(md_file, scale_factor, k_push, alp, n_runs, n_check, charge):
    """
    Do n_runs metadynamic runs with the specified parameters
    """
    directory_name = str(scale_factor)+'_'+str(k_push)+'_'+str(alp)
    os.mkdir(directory_name)
    shutil.copy('../'+md_file, directory_name)
    shutil.copy('../metadyn.inp', directory_name)
    os.chdir(directory_name)
    jobids, job_dictionary = submit_metadynamics_runs(md_file, k_push, alp, scale_factor,
                                                      range(n_runs))
    wait_for_jobs_to_finish_meta(jobids, job_dictionary, md_file, k_push, alp,
                                 scale_factor)
    #check the trajectory for reactions
    check_metadyn_runs(n_runs, n_check, charge, 'combined_dataframe.pkl')

    os.chdir('../')

if __name__ == "__main__":

    SCALE_FACTOR_LIST = [0.8]
    K_PUSH_LIST = [0.05]
    ALP_LIST = [0.3]

    REACTANT_IDX = 1
    N_CHECK = 10 #check every 10th saved structure in the calculated directory
    N_RUNS = 10 #do N_RUNS metadyn runs per parameter combination
    FILE_NAME = str(REACTANT_IDX) + '.xyz'
    SMILES = '[C:1](=[O:2])([C@:3]([C@:5]([O:8][O:11][H:12])([H:9])[H:10])([H:6])[H:7])[H:4]'
    MOL = Chem.MolFromSmiles(SMILES, sanitize=False)
    MOL = reorder_atoms_to_map(MOL)
    write_xyz_file(MOL, FILE_NAME)

    #Do md runs
    os.mkdir('md')
    shutil.copy(FILE_NAME, 'md')
    shutil.copy('md.inp', 'md')
    os.chdir('md')
    MD_IN_FILE = FILE_NAME
    for scale in SCALE_FACTOR_LIST:
        MD_OUT_FILE = calculate_md_relaxed_structure(MD_IN_FILE, scale, REACTANT_IDX)
        MD_IN_FILE = MD_OUT_FILE
    os.chdir('../')

    #Do metadynamics runs
    os.mkdir('metadynamics')
    os.chdir('metadynamics')

    for SCALE_FACTOR, K_PUSH, ALP in itertools.product(SCALE_FACTOR_LIST,
                                                       K_PUSH_LIST, ALP_LIST):
        MD_FILE = str(SCALE_FACTOR)+'_md.xyz'
        test_parameters(MD_FILE, SCALE_FACTOR, K_PUSH, ALP, N_RUNS, N_CHECK,
                        charge=0)
