import sys
import pandas as pd

job_df = pd.read_csv(sys.argv[1])

failed_runs = []

for idx in job_df.index:
    reactant_idx = job_df.loc[idx, 'reactant']
    r_idx = job_df.loc[idx, 'r_idx']
    letter = job_df.loc[idx, 'letter']
    pkl_file = str(reactant_idx)+'/'+str(reactant_idx)+'_'+str(r_idx)+'_'+letter+'.pkl'
    try:
        df = pd.read_pickle(pkl_file)
        e_hartree = df.loc[idx, "maxE"]
        job_df.loc[idx, 'maxE'] = e_hartree
        AB_hartree = df.loc[idx, "maxAB"]
        job_df.loc[idx, 'maxAB'] = AB_hartree
        BC_hartree = df.loc[idx, "maxBC"]
        job_df.loc[idx, 'maxBC'] = BC_hartree
        job_df.loc[idx, 'N_steps'] = df.loc[idx, 'N_steps']
    except FileNotFoundError:
        failed_runs.append(idx)

    print(idx)

print(job_df)
job_df.to_csv('updated_job_df.csv')

with open('failed_reactions.txt', 'w') as f:
    for item in failed_runs:
        f.write("%s\n" % item)
