# Bayesian Isotonic

This repository contains code for the methodology proposed in Boonstra, Owen, 
and Kang (2021). You can use the scripts in this repository to 
reproduce the simulations studies in the manuscript. See also `vignette.pdf` 
for a vignette describing the typical usage of the adaptive Bayesian priors 
on a simulated dataset

### Further details

In more detail, there are 11 files included in this repository (in addition to 
this README and `vignette.pdf`): one text file (`varying_data_run_sims.txt`), 
eight R scripts (ending in  `.R`), and two stan files (ending in `.stan`)

#### Text file

`varying_data_run_sims.txt` is the script for submitting parallel runs of
`varying_data_run_sims.R` (described below) to a cluster that is running
SLURM. The following command run in terminal will do this:

`sbatch varying_data_run_sims.txt`

The script assumes that you want all of the results to be put in your home 
directory (which you probably don't). Edit the script as needed  

#### `R` files

  - `functions.R` provides the R functions that implement the primary
methodological contribution of the paper. The main function is called
`bayesian_isotonic` and it can implement both the GAIPV and HSIPV priors from
the paper, depending on the Stan script provided as the value provided
to the argument `stan_path` (see Stan files section below)

  - `fixed_data_evaluation1.R` runs the first fixed-data evaluation in
Section 3.1 and creates Table 1 as LaTeX code

  - `fixed_data_evaluation2.R` runs the second fixed-data evaluation in
Section 3.2 and creates Figures 1 and 2

  - `varying_data_generate_params.R` constructs inputs for running the 
varying-data evaluation in Section 3.3. As described in the script's documentation 
and the language  below, these inputs can be overwritten by the user

  - `varying_data_run_sims.R` is the script to conduct the large-scale 
simulation study described in the manuscript. On a local machine, the user may 
choose a specific `array_id` (as described in this script's 
documentation) and run the code locally on their own machine. On a cluster 
running SLURM, the user can use this script to submit multiple jobs 
simultaneously (as described  in the description of `varying_data_run_sims.txt` above)

  - `varying_data_functions.R` provides the simulation functions. NOTE: 
this only contains the extra functions needed to simulate the data; the 
statistical methods are implemented in `functions.R`, described above

  - `process_varying_data_sims.R` gives the code to create the figures and tables in 
the manuscript and supplementary material reporting on the simulation study

  - `vignette.R` provides the R code for running the vignette yourself. 


#### Stan files

  - `iso_gamma.stan` implements the Gamma isotonic probability vector (GAIPV)
  prior 
  
  - `iso_horseshoe.stan` implements the Horseshoe isotonic probability vector
  (HSIPV) prior

#### Current Suggested Citation

Boonstra, Philip S.; Owen, Daniel R.; and Kang, Jian, "Shrinkage Priors for Isotonic Probability Vectors and Binary Data Modeling" (January 2020). The University of Michigan Department of Biostatistics Working Paper Series. Working Paper 128.
https://biostats.bepress.com/umichbiostat/paper128


