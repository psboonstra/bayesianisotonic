# Bayesian Isotonic

This repository contains code for the methodology proposed in Boonstra, Owen, and Kang (2021) You can use the scripts in this repository to 
reproduce the simulations studies in the manuscript. See also `vignette.pdf` 
for a vignette describing the typical usage of the adaptive Bayesian priors on
a simulated dataset

### Further details

In more detail, there are six files included in this repository (in addition to 
this README and `vignette.pdf`): one text file (<samp>run_abu_sims.txt</samp>) 
and five <samp>R</samp> scripts (ending in  <samp>.R</samp>). The simulation
studies reported in Boonstra and Barbaro were run using commit 22

#### Text file
<samp>run_abu_sims.txt</samp> is the script for submitting parallel runs of
<samp>run_aub_sims.R</samp> (described below) to a cluster that is running
SLURM. The following command run in terminal will do this:

<code> sbatch run_abu_sims.txt </code>

The script assumes that you want all of the results to be put in your home 
directory (which you probably don't). Edit the script as needed  

#### <samp>R</samp> files

<samp>vignette.R</samp> creates a single simulated dataset and walks through 
analyzing these data using the various adaptive priors. It can also be `knit` 
to create a copy of `vignette.pdf` (it will take a few minutes to knit)

<samp>0generate_params.R</samp> constructs inputs for running the simulation
study. As described in the script's documentation and the language below, these
inputs can be overwritten by the user

<samp>1run_iso_sims.R</samp> is the script to conduct the large-scale simulation study described in the manuscript. On a local machine, the user may choose a specific <samp>array_id</samp> (as described in this script's documentation) and run the code locally on their own machine. On a cluster running SLURM, the user can use this script to submit multiple jobs simultaneously (as described  in the description of <samp>run_abu_sims.txt</samp> above). 

<samp>2functions_simulation.R</samp> provides the simulation functions. NOTE: 
this only contains the extra functions needed to simulate the data; the methods
are contained in the next file below

<samp>3functions_methods.R</samp> provides the R functions that implement the functions themselves.

<samp>make_figures.R</samp> gives the code to create the figures and tables in 
the manuscript and supplementary material reporting on the simulation study. 

#### Current Suggested Citation

Boonstra, Philip S.; Owen, Daniel R.; and Kang, Jian, "Shrinkage Priors for Isotonic Probability Vectors and Binary Data Modeling" (January 2020). The University of Michigan Department of Biostatistics Working Paper Series. Working Paper 128.
https://biostats.bepress.com/umichbiostat/paper128


