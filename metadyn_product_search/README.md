# Meta-Dynamics search

## Usage 
Procedure to run meta-dynamics search for potential one-step reaction products:
* create .csv file containing the list of mapped reactant SMILES to be investigated
* change parameters CPUS, MEM, MAX_QUEUE, N_RUNS, S_FACTOR, TIME_PS, K_PUSH, ALP and SCRIPT in ```control_metadyn_runs.py```
  * CPUS: number of cpus each meta-dynamics job can use
  * MEM: amount of memory each meta-dynamics job can use
  * MAX_QUEUE: maximally allowed number of jobs to be submitted to ```Slurm``` at the same time
  * N_RUNS: number of meta-dynamics runs to be done for each reactant SMILES
  * S_FACTOR: the scaling parameter for the wall potential
  * K_PUSH: the k_push meta-dynamics parameter
  * ALP: parameter controlling the width of the Gaussian meta-dynamics potential 
  * SCRIPT: the path to ```reaction_box_no_slurm.py```
* change the path to xTB in ```reaction_box_no_slurm.py```


## External Programs:
The procedure relies on ```xTB``` for the quantum chemical calculations and ```Slurm``` for job submission
Also, a local version of ```xyz2mol``` (```xyz2mol_local.py```) is used to analyze the structures on the trajectory.
```
https://github.com/jensengroup/xyz2mol
```


Finally, the search is executed by:
```
./control_metadyn_runs.py reactant_smiles.csv
```

A folder named ```$S_FACTOR_$K_PUSH_$ALP``` is created, where the results are saved
Results are saved in  as dataframes in .pkl files: 1 for each run

The results for all runs can be combined using ```combine_runs.py```, which collects all recorded reactions and extracts all 1-step reactions, 
both are saved as .csv files.

## Test System
The folder ```test``` contains example output for cyclopropane created by running the following commands
```
mkdir test
cd test 
```
```reactant_smiles.csv``` is created containing the mapped SMILES for cyclopropane
```
../control_metadyn_runs.py reactant_smiles.csv
```
After the jobs have finished (change folder to fit your chosen parameter set and reactant index):
```
cd 0.6_0.03_0.7
python ../../combine_runs.py ../reactant_smiles.csv
cd 0
vi 0_1step.csv
```
which shows the recorded 1-step reactions for cyclopropane








