#!/groups/kemi/mharris/.conda/envs/rdkit_2020_09/bin/python

import os
import textwrap
import sys
import time

import pandas as pd


def qsub_prep(script_path, cpus, mem, smiles_idx, run_nr, smiles, s_factor,
              time_ps, k_push, alp):
    """
    write qsub file for SLURM subsmissin
    """
    pwd = os.getcwd()

    qsub_file = """\
    #!/bin/sh
    #SBATCH --job-name={3}_{4}
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task={1}
    #SBATCH --mem={2}
    #SBATCH --error={10}/{3}/run{4}.stderr
    #SBATCH --output={10}/{3}/run{4}.stdout
    #SBATCH --ntasks=1
    #SBATCH --time=8:00:00
    #SBATCH --partition=kemi1
    #SBATCH --no-requeue

    #mkdir /scratch/$SLURM_JOB_ID
    cd /scratch/$SLURM_JOB_ID

    #run python code
    ({0} {3} {4} '{5}' {6} {7} {8} {9})
    #cp output back

    cp run{4}/dataframe.pkl {10}/{3}/run{4}.pkl
    tar -zcvf structure_database.tar.gz run{4}/structure_database
    tar -zcvf run{4}.tar.gz run{4}

    cp run{4}.tar.gz {10}/{3}/run{4}.tar.gz
    cp structure_database.tar.gz {10}/{3}/run{4}_database.tar.gz

    rm {10}/{3}_run{4}_qsub.tmp

    #rm -r /scratch/$SLURM_JOB_ID

    """.format(script_path, cpus, mem, smiles_idx, run_nr, smiles, s_factor,
               time_ps, k_push, alp, pwd)

    with open(str(smiles_idx) + '_run'+str(run_nr) + "_qsub.tmp", "w") as qsub:
        qsub.write(textwrap.dedent(qsub_file))

    return str(smiles_idx) + '_run' + str(run_nr) + "_qsub.tmp"


def run_calculations(smiles_df, n_runs, script_path, s_factor, time_ps, k_push,
                     alp, max_queue, cpus, mem):
    """
    For each given smiles - submit n_runs metadynamics searches with the given
    parameters. Only submit new jobs when less than max_queue jobs in the queue
    """
    submitted_jobs = set()
    submitted_smiles = set()
    for idx in smiles_df.index:
        smiles = smiles_df.loc[idx, 'smiles']
        if idx not in submitted_smiles:
            os.mkdir(str(idx))
            submitted_smiles.add(idx)
        for run_nr in range(n_runs):
            qsub_name = qsub_prep(script_path, cpus, mem, idx, run_nr, smiles,
                                  s_factor, time_ps, k_push, alp)
            slurmid = os.popen("sbatch " + qsub_name).read()
            slurmid = int(slurmid.strip().split()[-1])

            submitted_jobs.add(slurmid)

            if len(submitted_jobs) >= max_queue:
                while True:
                    job_info = os.popen("squeue -u mharris").readlines()[1:]
                    current_jobs = {int(job.split()[0]) for job in job_info}

                    if len(current_jobs) >= max_queue:
                        time.sleep(30)
                    else:
                        finished_jobs = submitted_jobs - current_jobs
                        print("finished jobs: ", finished_jobs)
                        for job in finished_jobs:
                            submitted_jobs.remove(job)
                        break

    while True:
        job_info = os.popen("squeue -u mharris").readlines()[1:]
        current_jobs = {int(job.split()[0]) for job in job_info}
        if len(current_jobs) > 0:
            time.sleep(30)
        else:
            break


if __name__ == "__main__":
    SMILES_LIST = sys.argv[1]
    SMILES_DF = pd.read_csv(SMILES_LIST, index_col=0)
    print(SMILES_DF)
    CPUS = 1
    MEM = "3GB"
    MAX_QUEUE = 300
    N_RUNS = 100
    S_FACTOR = 0.6
    TIME_PS = 5
    K_PUSH = 0.03
    ALP = 0.7

    SCRIPT = '/groups/kemi/mharris/big_data_set/better_implementation/time_split_test/reaction_box_no_slurm.py'

    os.mkdir(str(S_FACTOR)+'_'+str(K_PUSH)+'_'+str(ALP))
    os.chdir(str(S_FACTOR)+'_'+str(K_PUSH)+'_'+str(ALP))

    run_calculations(SMILES_DF, N_RUNS, SCRIPT, S_FACTOR, TIME_PS, K_PUSH,
                     ALP, MAX_QUEUE, CPUS, MEM)

    #collect dataframes
    #dict_rows = []
    #for ridx in JOB_NAMES:
    #    df = pd.read_pickle(str(ridx)+'.pkl')
    #    rows_dict = df.to_dict()
    #    print(rows_dict)
    #    dict_rows.append(rows_dict)

    #UPDATED_DF = pd.DataFrame(dict_rows)
    #print(UPDATED_DF)
    #UPDATED_DF.to_pickle('updated_systematic_rmsd_v2_run1.pkl')
