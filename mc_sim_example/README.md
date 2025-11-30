# Monte Carlo simulation sequence designs
## Overview and installation
An example script for designing sequences from the Potts model coefficients for a given protein family is found in the `potts_mc_sim.py` script. A conda environment for running the script is specified in the `mc_environment.yml` file. The conda environment can be installed by:
```
conda env create -f mc_environment.yml
conda activate mc_env
```

## Run parameters
### Input parameters
Input parameters for a design run are specified in the script at the top of the script. Input parameters to set are: 
 * `path`: Directory path for running simulation and saving outputs.
 * `msaFile`: MSA file for protein family used for fitting the Potts model.
 * `hFile`: Single-site (h) coefficients from the fit of the Potts model for the protein family from the scripts in the ../plmDCA_asymmetric_v2_fit directory.
 * `jFile`: Coupling (J) coefficients from the fit of the Potts model for the protein family from the scripts in the ../plmDCA_asymmetric_v2_fit directory.
 * `jMatName`: Name for the gauge-transformed J-coefficient matrix in the Potts model fit matLab script. Should be "Jtrans" for script(s) used in the ../plmDCA_asymmetric_v2_fit directory.

 ## Simulation hyperparamters
 Simulation hyperparameters for a design run are specified in the script at the top of the script. Hyperparameter to set are:
  * `seed`: Value for seeding random.
  * `nRounds`: Number of rounds for the Monte Carlo simulation.
  * `nSeqMSAs`: Number of sequences to design in independent trajectories.
  * `betaSet`: Initial value for beta inverse temperature parameter. 
  * `betaInc`: Value to increase the beta inverse temperature parameter (i.e. lower the simulation "temperature") on all increase steps.
  * `betaStepInc`: The number of simulation rounds at each value of beta. As specified, parameter is set to increase for a number of rounds such that the simulation contains 20 increases in beta. 
  * `hWeight`: Weight for H(Seq) in the sequence energy function E(seq) = hWeight * H(seq) + jWeight * J(seq). A value of 1 is used for all simulaitions unless otherwise specified.
  * `jWeight`: Weight for J(Seq) in the sequence energy function E(seq) = hWeight * H(seq) + jWeight * J(seq). A value of 1 is used for all simulaitions unless otherwise specified.
  * `numCores`: Number of workers to use for multiprocessing. 

# Running simulations
A Monte Carlo simulation for designing sequences from the Potts model coefficients for a given protein family can be run by:

``` python potts_mc_sim.py```

The script with run the simulation with the parameters specified above. The simulation will result in an ensemble of sequences designed by the specified energy function. The final designed sequences at the final round of the simulation will be found in the `final_seqs.txt` file. A sequence logo for the final sequences can be generated using [Weblogo](https://weblogo.berkeley.edu/logo.cgi). 

Figures overview the simulation run will be saved in the defined path directory, including:
 * `MC_Etot.png`: E(seq) energies for all sequence design trajectories throughout the simulation.
 * `MC_htot.png`: H(seq) energies for all sequence design trajectories throughout the simulation.
 * `MC_Jtot.png`: J(seq) energies for all sequence design trajectories throughout the simulation.