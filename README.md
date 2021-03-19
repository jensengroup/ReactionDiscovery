# ReactionDiscovery

## Systematic Search


## Meta-dynamics Search
Searches the product space based on meta-dynamics runs

### Usage
To run meta-dynamics search in order to locate potential one-step products:
* change the paramater SMILES in ```reaction_box.py``` to the mapped SMILES of your reactant
* change the parameters SCALE_FACTOR_LIST, K_PUSH_LIST and ALP_LIST to the meta-dynamics parameters you want to test
* change N_RUNS to the number of times you want each parameter set to be run
* ```md.inp``` and ```metadyn.inp``` should be present in the working directory
  * if different parameters relating to the MD run are wanted; change these.


### External programs
The procedure relies on ```xTB``` for the quantum chemical calculations. The are submitted through the submit scripts localted in the ```submit_scripts``` folder
Also, a local version of ```xyz2mol``` (```xyz2mol_local.py```) is used to analyze the structures on the trajectory.
```
https://github.com/jensengroup/xyz2mol
```

Finally, the search is executed by
```
./reaction_box.py
```

Results will we found in ```metadynamics``` containing one folder for each parameter setting. In each og these folders a file called ```combined_dataframe.pkl``` will contain information about the reactions located during the search.
