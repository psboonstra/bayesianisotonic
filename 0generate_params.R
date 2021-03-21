#DESCRIPTION: creates a list of lists called 'arglist'. The user should not generally change anything in this script. If it is desired to overwrite
#any of the default values, create that variable in 'runMe.R' (before sourcing GenParams.R on line 65 of that script). The default value below will 
#then be ignored. 

if(!"n_mc_warmup"%in%ls()){n_mc_warmup = 2.5e3;}
if(!"n_mc_samps"%in%ls()){n_mc_samps = 2.5e3;}
if(!"mc_chains"%in%ls()){mc_chains = 2;}
if(!"mc_thin"%in%ls()){mc_thin = 1;} 
if(!"mc_stepsize"%in%ls()){mc_stepsize = 0.1;}
if(!"mc_adapt_delta"%in%ls()){mc_adapt_delta = 0.99;}
if(!"mc_max_treedepth"%in%ls()){mc_max_treedepth = 18;}
if(!"ntries_per_iter"%in%ls()){ntries_per_iter = 1;}
if(!"include_nonbayes"%in%ls()){include_nonbayes = TRUE;}
if(!"ga_lower_bound"%in%ls()){ga_lower_bound = .Machine$double.eps;}
if(!"return_summarized_models"%in%ls()){return_summarized_models = TRUE;}
if(!"nonbayes_added_weight"%in%ls()){nonbayes_added_weight = c(0, 1/64)}

if(!"dynamic_run"%in%ls()){dynamic_run = F;}

if(!"hs_stan_filenames"%in%ls()) {
  hs_stan_filenames = c(horseshoe = "iso_horseshoe.stan");
}
if(!"ga_stan_filenames"%in%ls()) {
  ga_stan_filenames = c(gamma = "iso_gamma.stan");
}

if(!"n_training_seq"%in%ls()){
  n_training_seq = c(80, 320);
}

if(!"predictor_dist_list"%in%ls()) {
  predictor_dist_list = list(
    #Dist 1: uniform discrete 
    function(n) {sample(x = seq(0.05, 0.95, length = 10), 
                        size = n, 
                        replace = TRUE)}
  )
}

if(!"true_prob_curve_list"%in%ls()) {
  true_prob_curve_list = list(
    # Curve 1: simple linear increase
    function(x) {pmin(1, pmax(0, x));},
    # Curve 2: Two smooth, sharp increases:
    function(x) {
      0.35 * plogis(x, location = 0.25, scale = 0.015) +
        0.65 * plogis(x, location = 0.75, scale = 0.015);
    }
  )
}

if(!"n_sim"%in%ls()) {
  n_sim = 10;
}

all_varying = 
  crossing(n_training = n_training_seq, 
           predictor_dist = 1:length(predictor_dist_list),
           true_prob_curve = 1:length(true_prob_curve_list)) 

all_varying = cbind(all_varying,
                    scenario = (1:nrow(all_varying)));

arglist = list();
for(i in 1:nrow(all_varying)) {
  

  foo = list(array_id = array_id,
             scenario_id = all_varying[i, "scenario"],
             n_sim = n_sim,
             n_training = all_varying[i, "n_training"],
             n_validation = 1e3,
             true_prob_curve = true_prob_curve_list[[all_varying[i, "true_prob_curve"]]],
             true_prob_curve_id = all_varying[i, "true_prob_curve"],
             predictor_dist = predictor_dist_list[[all_varying[i, "predictor_dist"]]],
             predictor_dist_id = all_varying[i, "predictor_dist"],
             hs_stan_filenames = hs_stan_filenames,
             ga_stan_filenames = ga_stan_filenames,
             ga_lower_bound = ga_lower_bound,
             n_mc_warmup = n_mc_warmup, 
             n_mc_samps = n_mc_samps, 
             mc_chains = mc_chains, 
             mc_thin = mc_thin,
             mc_stepsize = mc_stepsize,
             mc_adapt_delta = mc_adapt_delta,
             mc_max_treedepth = mc_max_treedepth,
             ntries_per_iter = ntries_per_iter,
             random_seed = NULL,
             data_seeds = NULL,
             stan_seeds = NULL,
             include_nonbayes = include_nonbayes,
             return_summarized_models = return_summarized_models,
             nonbayes_added_weight = nonbayes_added_weight,
             dynamic_run = dynamic_run);
  arglist = c(arglist, list(foo));
}
