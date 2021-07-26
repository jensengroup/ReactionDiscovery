#!/groups/kemi/mharris/.conda/envs/rdkit_2020_09/bin/python

import shutil
import os
import tarfile
import numpy as np
import pandas as pd


def check_path_interpolation(directory):
    os.chdir(directory)
    tfs = [f for f in os.listdir(os.curdir) if \
            f.endswith("log_xtb.tar.gz")]

    for file_name in tfs:
        tf = tarfile.open(file_name, mode='r')
        tf.extractall()
        tf.close()

    files = [f for f in os.listdir(os.curdir) if \
             f.endswith("xtbout")]
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    max_energy = None
    for file_name in files:
        with open(file_name, 'r') as _file:
            line = _file.readline()
            while line:
                if 'TOTAL ENERGY' in line:
                    energy_au = np.float(line.split()[3])
                line = _file.readline()
        if not max_energy:
            max_energy = energy_au
            ts_guess = file_name[:-6]+'xyz'
        if energy_au > max_energy:
            max_energy = energy_au
            ts_guess = file_name[:-6]+'xyz'
        #os.remove(file_name)
    print(ts_guess, max_energy)

    os.chdir('../')
    return ts_guess, max_energy



def find_ts_guess(directory, rxn_index, letter):
    """
    when sp calculations ar efinished: find the structure with maximum xtb
    energy
    """
    ts_guess_paths = []
    ts_guess_energies = []
    os.chdir(directory)
    high_temperature = False
    if os.path.exists('ht'):
        os.chdir('ht')
        high_temperature = True

    paths = [d for d in os.listdir(os.curdir) if d.startswith('path') and os.path.isdir(d)]
    print(paths)
    for path in paths:
        ts_guess, max_energy = check_path_interpolation(path)
        ts_guess_paths.append(path+'/'+ts_guess)
        ts_guess_energies.append(max_energy)


    ts_guess = ts_guess_paths[ts_guess_energies.index(max(ts_guess_energies))]
    if high_temperature:
        os.chdir('../')
        ts_guess = 'ht/'+ts_guess

    shutil.copy(ts_guess, '../../dft_ts/'+rxn_index+'_'+letter+'.xyz')

    os.chdir('../')


if __name__ == "__main__":
    df = pd.read_pickle('updated_dataframe5.pkl')
    os.chdir('push_pull5')
    reactant_xtb_energy = -21.569104 #Hartree
    letter = 'e'
    for rxn in df.index:
        barrier_estimate = (df.loc[rxn, 'xTB max [Hartree]']-reactant_xtb_energy)*627.509
        if barrier_estimate < 50:
            rxn_index = (df.loc[rxn, 'reac_num'])[5:]
            print(rxn, rxn_index, barrier_estimate)
            find_ts_guess(str(rxn), rxn_index, letter)
