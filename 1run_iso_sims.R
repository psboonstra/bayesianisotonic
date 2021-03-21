# DESCRIPTION: Running this script without modification produces an object 
# called 'arglist', which is a list of 4 lists, with each of the 4 lists 
# corresponding to a unique scenario, defined as a combinations of 2 sample sizes
# sizes and 2 probability curves. Then, one of these 
# lists is extracted based upon the value of array_id, and the simulator function 
# called simulator is called on this scenario. 

# This script expects that the following files are in your current working directory:
# iso_horseshoe.stan, 
# iso_gamma.stan, 
# generate_params.R
# functions.R
# You can point to a different directory for the stan files by changing the object 'stan_file_path'. 


library(tidyverse);library(Iso);library(cir);library(binom);

#library(bisoreg);
#Flag for whether this is running on a local machine or on a cluster running SLURM
my_computer = F;

#The simulation study in Boonstra and Barbaro was conducted in two batches of 400 jobs followed by another 400 (because the cluster only accepts 
#jobs in batches of size up to 1000). 'which_run' indicates whether this is the first batch (= 1 ) or the second batch ( = 2)
which_run = 2;

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

# 'jobs_per_scenario' is the number of parallel independent jobs to send for 
# each scenario, and 'n_sim' is the number of independent datasets per job
# to generate. The constraints are as follows: 
# (i) jobs_per_scenario * length(arglist) < 1000 (because only 1000 jobs
# can be submitted at once)
# (ii) sum_{i=1}^{jobs_per_scenario} n_sim_i = total desired sims per scenario

jobs_per_scenario = ifelse(my_computer, 1, 100);
# Should the sim_id labels be randomly permuted across array_ids?
permute_array_ids = ifelse(my_computer, F, T);

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
# sourcing generate_params.R 
source("generate_params.R");

rm(list=setdiff(ls(all=T),c("arglist","array_id","my_computer",
                            "permute_array_ids","jobs_per_scenario",
                            "array_id_offset")));
source("functions_simulation.R");
source("functions_methods.R");

#This is the step that permutes the duplicated job ids (the first are left alone)
#If FALSE, then all jobs from the same scenario will occur in contiguous blocks. 
if(permute_array_ids) {
  permute_array_ids = seq(1, jobs_per_scenario * length(arglist), by = jobs_per_scenario)
  set.seed(2);
  permute_array_ids = 
    c(permute_array_ids, 
      setdiff(sample(jobs_per_scenario * length(arglist)), permute_array_ids));
  
} else {
  permute_array_ids = 1:(jobs_per_scenario * length(arglist));
}

curr_args = arglist[[ceiling(permute_array_ids[array_id]/jobs_per_scenario)]];

curr_args[["random_seed"]] = array_id;


# Phil delete below this line! ----
#curr_args$n_sim = 8;
#curr_args$random_seed = 200;
#curr_args$data_seeds = 794080207
#curr_args$stan_seeds = 2046114256
#curr_args$n_training = 100;
#curr_args$ga_stan_filenames = NA;
#  c(gamma3 = "iso_gamma.stan");
#curr_args$ga_lower_bound = 
#  c(gamma3 = .Machine$double.eps^(1.25));
#curr_args$include_nonbayes = FALSE;
#curr_args$do_these_cut_strategies = NA;
#curr_args$true_prob_curve =  
#  function(x, eps = 0.01) {eps + (1 - 2 * eps) * (x > 0.25);}
#curr_args$include_nonbayes = FALSE;
#curr_args$ga_stan_filenames = NA;
#curr_args$do_these_cut_strategies = c("pava");
#curr_args$n_mc_samps = 250;
#curr_args$n_mc_warmup = 250;
#curr_args$random_seed = 850030880;
#stop();
# Phil delete above this line! ----

assign(paste0("job",array_id_offset + array_id),
       do.call("simulator",args = curr_args));

if(0) {
  foo <- 
    get(paste0("job",array_id_offset + array_id))$summarized_performance;
  
  foo %>% 
    group_by(priors, cut_strategies) %>%
    summarize(rmse = mean(rmse), 
              bias = mean(bias), 
              kl_div = mean(kl_div),
              coverage_50 = mean(coverage_50), 
              loglik = mean(loglik)) %>%
    as.data.frame()
  
  
  foo2 <- 
    get(paste0("job",array_id_offset + array_id))$summarized_models
  
  
  foo2 %>%
    filter(x %in% x[c(1,20, 50,80, 101)]) %>%
    group_by(priors, cut_strategies, x) %>%
    summarize(true_prob = mean(true_prob), 
              mean_fitted_prob = mean(fitted_prob), 
              median_fitted_prob = median(fitted_prob));
  
  ggplot() +
    geom_line(data = foo2,
              aes(x = x, 
                  y = fitted_prob, 
                  color = priors, 
                  group = interaction(priors, array_id, sim_id)),
              alpha = 0.9, 
              size = 0.4) +
    geom_line(data =
                filter(foo2,
                       predictor_dist_id == 1,
                       #priors == dplyr::first(priors),
                       #cut_strategies == dplyr::first(cut_strategies),
                       sim_id == 1) %>%
                arrange(true_prob_curve_id, array_id) %>%
                group_by(true_prob_curve_id, n_training) %>%
                filter(array_id == first(array_id)),
              aes(x = x, 
                  y = true_prob), 
              color = "black",
              alpha = 1,
              size = 0.6) +
    facet_grid(n_training ~ cut_strategies, scales = "free_y") + 
    coord_cartesian(xlim = c(-2, 2), ylim = c(0,1)) + 
    scale_y_continuous(name = "Fitted / true probability") + 
    scale_x_continuous(name = "X", labels = NULL) + 
    scale_color_brewer(name = "Prior",
                       palette = "Set2") + 
    theme(legend.position = "top", 
          legend.direction = "horizontal", 
          text = element_text(size = 16));
  
}

do.call("save",list(paste0("job",array_id_offset + array_id),
                    file = paste0("out/job",array_id_offset + array_id,".RData"),
                    precheck = FALSE));

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