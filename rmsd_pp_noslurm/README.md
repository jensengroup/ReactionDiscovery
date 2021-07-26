# RMSD-PP procedure

If products are from the meta-dynamics product search, reactant and product structures can be extracted from the meta-dynamics
search. This option is set by calling --use-structures when creating the submission file in ```rmsd_pp_no_slurm.py```.
Otherwise structures are embedded using ```RDKit```. 

## Usage

A .csv file containing the reactions to get barrier estimates mus be created: these reactions can stem from either 
meta-dynamics search, systematic search or something else entirely.
The .csv file must contain one row for each reaction with at least two columns: 
* ```reactant_smiles_am``` containing the mapped SMILES of the reactant
* ```product_smiles_am``` containing the mapped SMILES of the product
Then change parameters CPUS, MEM, MAX_QUEUE and SCRIPT in ```control_rmsd_paths.py``` 
* CPUS: number of cpus each meta-dynamics job can use
* MEM: amount of memory each meta-dynamics job can use
* MAX_QUEUE: maximally allowed number of jobs to be submitted to ```Slurm``` at the same time
* SCRIPT: the path to ```rmsd_pp_no_slurm.py```

### for meta-dyanmics products
```create_job_names.py``` can be used to make the reaction.csv file after meta-dynamics runs. 
By default 5 entries per reaction is created.
RXN_INDEX_LIST must be changed to the list of reactant indexes used in the meta-dynamics search


The RMSD-PP paths are then computed by running
```
./control_rmsd_paths.py reactions.csv
```

## External Programs:
The procedure relies on ```xTB``` for the quantum chemical calculations and ```Slurm``` for job submission
Also, a local version of ```xyz2mol``` (```xyz2mol_local.py```) is used to analyze the structures on the trajectory.
```
https://github.com/jensengroup/xyz2mol
```

Results for each run is saved in a .pkl file named according to the job name.

The results for all runs can be collected by 
```
./collect_dataframes.py reactions.csv
```
