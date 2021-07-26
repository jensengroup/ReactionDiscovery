# Finding and testing DFT transition states (TSs)
This procedures optimizes a TS based on the specified DFT functional and basis set and guess structures from e.g. the RMSD-PP procedure.

## Usage
As input is needed the results from the RMSD-PP procedure: both the paths containing the TS guesses, and the updated .csv file containing barrier estimates for the reactions. 

Change the parameters cpus, mem, max_queue and script in ```control_dft_ts.py``` to fit your requirements
* cpus: number of cpus each job should use
* mem: amount of memory each job should use
* max_queue: maximally allowed number of submitted jobs at a time
* script_path: path to ```run_dft_ts.py```

remember to change the paths to xTB and Gaussian.


## External Programs:
The procedure relies on ```Gaussian``` for the quantum chemical calculations and ```Slurm``` for job submission
Also, a local version of ```xyz2mol``` (```xyz2mol_local.py```) is used to analyze the structures on the trajectory.
```
https://github.com/jensengroup/xyz2mol
```

The calculations are done by calling
```
./control_dft_ts.py rmsd_output.csv
```
