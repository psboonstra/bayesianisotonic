# DESCRIPTION: Running this script without modification produces an object 
# called 'arglist', which is a list of 4 lists, with each of the 4 lists 
# corresponding to a unique scenario, defined as a combinations of 2 sample sizes
# sizes and 2 probability curves. Then, one of these 
# lists is extracted based upon the value of array_id, and the simulator function 
# called simulator is called on this scenario. 
# 
# This script expects that the following files are in your current working directory:
# iso_horseshoe.stan, 
# iso_gamma.stan, 
# varying_data_generate_params.R
# varying_data_functions.R
# functions.R

library(tidyverse);library(Iso);library(cir);library(binom);

#Flag for whether this is running on a local machine or on a cluster running SLURM
my_computer = F;

if(my_computer) {
  library(rstan);
  #Choose from between 1-400 if generate_params.R is used as-is
  array_id = 1;
  rstan_options(auto_write = TRUE);
} else {
  library(rstan);
  array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));
  rstan_options(auto_write = FALSE);
}

#Helpful to print warnings when they occur for debugging
options(warn = 1);

#Recommended options from rstan:
options(mc.cores = parallel::detectCores());



# The varying-data evaluation in Section 3.3 was conducted in two batches of 
# 400 jobs each. For i=1,...,400, job i in each batch is seeded identically 
# and so samples identical datasets. The batches differ in the methods that they
# implement: the first batch implements the fast methods (HS and GA1) and the
# second batch implements the slow methods (GA2, GA3, GA4). Each job corresponds
# to one of the four true data generating scenarios described in Section 3.3, thus
# there are 400/4 = 100 jobs per scenario ('jobs_per_scenario = 100'). 
# Each job conducts two independent iterations ('nsim=2'), hence 200 iterations
# per scenario. 

jobs_per_scenario = ifelse(my_computer, 1, 100);
# Should the sim_id labels be randomly permuted across array_ids?
permute_array_ids = ifelse(my_computer, F, T);

# 'which_run' indicates whether this is the first batch (= 1 ) or the second batch ( = 2)
which_run = 2;

if(which_run == 1) {
  
  # Request 59 minutes running time for this batch (using 2 nodes)
  
  n_sim = ifelse(my_computer, 1, 2);
  
  # 'array_id_offset' is added to the label of the saved object. Useful when a 
  # new batch of jobs is run and you want to continue the labeling scheme. 
  array_id_offset = 0;
  
  # This is the main simulation study;
  hs_stan_filenames = c(horseshoe = "iso_horseshoe.stan");
  ga_stan_filenames = c(gamma1 = "iso_gamma.stan");
  ga_lower_bound = c(gamma1 = .Machine$double.eps);
  
} else if(which_run == 2) {#second batch
  
  # Request 2 days running time for this batch (using 1 node)
  
  n_sim = ifelse(my_computer, 1, 2);
  
  # number of jobs from first batch: jobs_per_scenario * length(arglist)
  array_id_offset = (jobs_per_scenario * 4);
  
  # This batch studies the effect of pushing the gamma lower bound closer to zero
  hs_stan_filenames = NA;
  ga_stan_filenames = c(gamma2 = "iso_gamma.stan", 
                        gamma3 = "iso_gamma.stan", 
                        gamma4 = "iso_gamma.stan");
  ga_lower_bound = c(gamma2 = 10^(log10(.Machine$double.eps) - 1), 
                     gamma3 = 10^(log10(.Machine$double.eps) - 2), 
                     gamma4 = 0);
  include_nonbayes = FALSE;
  
} 

# Before calling the next line, you should specify any parameters that 
# you wish to change from their default values. Any specified values will
# take precedence over the default values that will otherwise be et by
# sourcing varying_data_generate_params.R 
source("varying_data_generate_params.R");

rm(list=setdiff(ls(all=T),c("arglist","array_id","my_computer",
                            "permute_array_ids","jobs_per_scenario",
                            "array_id_offset")));
source("varying_data_functions.R");
source("functions.R");

# This is the step that permutes the duplicated job ids (the first are left alone)
# Do this if you plan to check the results along the way and want to ensure
# a good representation of all scenarios
# If FALSE, then all jobs from the same scenario will occur in contiguous blocks. 
if(permute_array_ids) {
  permute_array_ids = seq(1, jobs_per_scenario * length(arglist), by = jobs_per_scenario)
  set.seed(2);
  permute_array_ids = 
    c(permute_array_ids, 
      setdiff(sample(jobs_per_scenario * length(arglist)), permute_array_ids));
  
} else {
  permute_array_ids = 1:(jobs_per_scenario * length(arglist));
}

# It's wasteful, but we only need the single id for this job
curr_args = arglist[[ceiling(permute_array_ids[array_id]/jobs_per_scenario)]];
# Ensure that the datasets are identical within array_ids and different 
# between array_ids
curr_args[["random_seed"]] = array_id;

# This is the actual call to the simulator function
assign(paste0("job",array_id_offset + array_id),
       do.call("simulator",args = curr_args));

# Save the entire workspace in case you want to dig in
do.call("save",list(paste0("job",array_id_offset + array_id),
                    file = paste0("out/job",array_id_offset + array_id,".RData"),
                    precheck = FALSE));
# Also just save the performance as a csv 
write_csv(get(paste0("job",array_id_offset + array_id))$summarized_performance,
          path = paste0("out/job",array_id_offset + array_id,"_performance.csv"),
          append = FALSE);

write_csv(get(paste0("job",array_id_offset + array_id))$summarized_bayesian_performance,
          path = paste0("out/job",array_id_offset + array_id,"_bayesian_performance.csv"),
          append = FALSE);

if(curr_args[["return_summarized_models"]]) {
  write_csv(get(paste0("job",array_id_offset + array_id))$summarized_models,
            path = paste0("out/job",array_id_offset + array_id,"_models.csv"),
            append = FALSE);
}


if(0) {
  
  for(j in 1:length(true_prob_curve_list)) {
    for(i in 1:length(predictor_dist_list)) {
      cat(mean(true_prob_curve_list[[j]](predictor_dist_list[[i]](1e6))), "\n")
    }
  }
  
}