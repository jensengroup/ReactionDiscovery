# Finding and testing xTB transition states (TSs)
This procedures optimizes a TS based on a RMSD-PP guess structure using TS optimization algorithm in Gaussian, computing forces from GFN2-xTB

## Usage
As input is needed the results from the RMSD-PP procedure: both the paths containing the TS guesses, and the updated .csv file contining barrier estimates for the reactions. 

Change the parameters cpus, mem, max_queue and script in ```control_xtb_ts.py``` to fit your requirements
* cpus: number of cpus each job should use
* mem: amount of memory each job should use
* max_queue: maximally allowed number of submitted jobs at a time
* script: path to ```run_xtb_ts.py```

remember to change the paths to xTB and Gaussian.


## External Programs:
The procedure relies on ```xTB``` as well as ```Gaussian``` for the quantum chemical calculations and ```Slurm``` for job submission
Also, a local version of ```xyz2mol``` (```xyz2mol_local.py```) is used to analyze the structures on the trajectory.
```
https://github.com/jensengroup/xyz2mol
```

The calculations are done by calling
```
./control_xtb_ts.py rmsd_output.csv
```
