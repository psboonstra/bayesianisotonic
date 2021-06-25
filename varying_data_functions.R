
rgamma_trunc = function(n, shape, rate, lb = 0) {
  qgamma(runif(n) * pgamma(lb, shape, rate, lower = FALSE), shape, rate, lower = FALSE);
}


validate_bayesian_isotonic = function(x_validation_category, 
                                      true_prob_validation,
                                      draws_xi) {
  
  matrix_fitted_probs = 
    draws_xi[,x_validation_category,drop = FALSE];
  if(any(is.nan(matrix_fitted_probs))) {
    matrix_fitted_probs = 
      matrix_fitted_probs[which(rowSums(is.nan(matrix_fitted_probs)) == 0), , drop = F]
  }
  tiny_positive = .Machine$double.eps;
  matrix_fitted_probs = pmax(pmin(matrix_fitted_probs, 1 - tiny_positive), tiny_positive);
  
  matrix_true_probs = 
    matrix(rep(true_prob_validation, each = nrow(matrix_fitted_probs)), 
           nrow = nrow(matrix_fitted_probs));
  
  matrix_true_probs = pmax(pmin(matrix_true_probs, 1 - tiny_positive), tiny_positive);
  
  matrix_bias_probs = 
    matrix_true_probs - matrix_fitted_probs;
  
  curr_rmse = 
    sqrt(mean((matrix_bias_probs)^2));
  
  curr_bias = 
    mean(matrix_bias_probs);
  
  curr_loglik = 
    mean((matrix_true_probs * log(matrix_fitted_probs)) + 
           ((1 - matrix_true_probs) * log(1 - matrix_fitted_probs)) - 
           tiny_positive * log(tiny_positive));
  
  curr_kl_div = 
    mean((matrix_true_probs * log(matrix_true_probs)) + 
           ((1 - matrix_true_probs) * log(1 - matrix_true_probs)) -
           tiny_positive * log(tiny_positive)) - 
    curr_loglik;
  
  curr_coverage = 
    mean((apply(matrix_bias_probs, 2, quantile, p = 0.25) < 0) &
           ((apply(matrix_bias_probs, 2, quantile, p = 0.75) > 0)))
  
  summary_fitted_isotonic = 
    tibble(rmse = curr_rmse, 
           bias = curr_bias,
           loglik = curr_loglik,
           kl_div = curr_kl_div,
           coverage = curr_coverage);
  
  returned_predictions =
    tibble(x_cat = x_validation_category, 
           true_prob = true_prob_validation,
           fitted_prob = apply(matrix_fitted_probs, 2, median),
           lower50 = apply(matrix_fitted_probs, 2, quantile, p = 0.25),
           upper50 = apply(matrix_fitted_probs, 2, quantile, p = 0.75))
  
  
  list(summary_fitted_isotonic = summary_fitted_isotonic,
       returned_predictions = returned_predictions)
}


validate_isotonic = function(x_validation, 
                             true_prob_validation,
                             fitted_isotonic) {
  
  tiny_positive = .Machine$double.eps;
  
  fitted_isotonic <- 
    tibble(x_validation = x_validation, 
           true_prob_validation = true_prob_validation) %>%
    group_by(x_validation) %>%
    summarize(true_prob_validation = mean(true_prob_validation),
              number_validate = n()) %>%
    full_join(fitted_isotonic, 
              by = c("x_validation" = "x")) %>%
    arrange(x_validation) %>%
    fill(y, lower50conf, upper50conf, 
         .direction = "downup") %>%
    mutate(fitted_prob_validation =
             pmax(pmin(y, 1 - tiny_positive), tiny_positive),
           true_prob_validation =
             pmax(pmin(true_prob_validation, 1 - tiny_positive), tiny_positive), 
           bias_prob = 
             true_prob_validation - fitted_prob_validation,
           loglik = 
             (true_prob_validation * log(fitted_prob_validation)) + 
             ((1 - true_prob_validation) * log(1 - fitted_prob_validation)) - 
             tiny_positive * log(tiny_positive), 
           kl_div = 
             (true_prob_validation * log(true_prob_validation)) + 
             ((1 - true_prob_validation) * log(1 - true_prob_validation)) -
             tiny_positive * log(tiny_positive) - 
             loglik, 
           coverage = 
             (lower50conf <= true_prob_validation) * 
             (upper50conf >= true_prob_validation));
  
  
  summary_fitted_isotonic <- 
    fitted_isotonic %>%
    summarize(rmse = 
                sqrt(weighted.mean((bias_prob)^2, 
                                   w = number_validate)),
              bias = 
                weighted.mean(bias_prob, 
                              w = number_validate),
              loglik = 
                weighted.mean(loglik, 
                              w = number_validate),
              kl_div = 
                weighted.mean(kl_div, 
                              w = number_validate), 
              coverage = 
                weighted.mean(coverage, 
                              w = number_validate));
  
  list(fitted_isotonic = fitted_isotonic,
       summary_fitted_isotonic = summary_fitted_isotonic);
}




simulator = function(array_id = 1,
                     scenario_id = 1,
                     n_sim = 10,
                     n_training = 250,
                     n_validation = 2000,
                     true_prob_curve = function(x) {0.25 + 0.50 * (x > 0)}, 
                     true_prob_curve_id = 1,
                     predictor_dist = function(n) {rnorm(n = n);},
                     predictor_dist_id = 1,
                     hs_stan_filenames = 
                       c(horseshoe = "iso_horseshoe.stan"),
                     ga_stan_filenames = 
                       c(gamma = "iso_gamma.stan"),
                     ga_lower_bound = 
                       c(gamma = .Machine$double.eps),
                     n_mc_warmup = 2.5e3, 
                     n_mc_samps = 5e3, 
                     mc_chains = 1, 
                     mc_thin = 1, 
                     mc_stepsize = 0.1, 
                     mc_adapt_delta = 0.99,
                     mc_max_treedepth = 15,
                     ntries_per_iter = 1,
                     random_seed = sample.int(.Machine$integer.max,1),
                     data_seeds = NULL,#If non-null, should have length equal to n_sim
                     stan_seeds = NULL,#If non-null, should have length equal to n_sim
                     include_nonbayes = TRUE,
                     return_summarized_models = TRUE,
                     nonbayes_added_weight = c(0, 1/64), 
                     dynamic_run = T) 
{ 
  require(tidyverse);require(rstan);require(Iso);require(cir);require(binom);
  
  begin_all = Sys.time();
  set.seed(random_seed);
  if(!length(data_seeds)) {
    data_seeds = sample.int(.Machine$integer.max,n_sim);
  } else {
    if(length(data_seeds) != n_sim) {
      stop("'data_seeds' was provided but must have length equal to 'n_sim'");
    }  
  }
  if(!length(stan_seeds)) {
    stan_seeds = sample.int(.Machine$integer.max,n_sim);
  } else {
    if(length(stan_seeds) != n_sim) {
      stop("'stan_seeds' was provided but must have length equal to 'n_sim'");
    }  
  }
  
  # At least one of hs_stan_filenames or ga_stan_filenames must have a valid
  # entry, meaning a named, valid path. If both are NA, the simulator will
  # be a dry run with no actual methods fit
  
  if(is.na(hs_stan_filenames)[1] &&
     is.na(ga_stan_filenames)[1]) {
    do_these_priors = NULL;
    fit_methods = include_nonbayes;
  } else {
    fit_methods = T;
    if(any(is.na(hs_stan_filenames)) || !length(hs_stan_filenames)) {
      hs_stan_filenames = NULL;
    } else if(is.null(names(hs_stan_filenames))) {
      names(hs_stan_filenames) = 
        paste0("horseshoe",1:length(hs_stan_filenames));
    }
    
    # If ga_stan_filenames is missing, then skip it; 
    if(any(is.na(ga_stan_filenames)) || !length(ga_stan_filenames)) {
      ga_stan_filenames = NULL;
    } else if(is.null(names(ga_stan_filenames))) {
      names(ga_stan_filenames) = 
        paste0("gamma",1:length(ga_stan_filenames));
    }
    
    # if ga_stan_filenames was provided, it must have accompanying lower bounds
    if(!is.null(ga_stan_filenames))  {
      if(length(ga_lower_bound) == 1) {
        # recycle lower bound as needed
        ga_lower_bound = rep(ga_lower_bound, length = length(ga_stan_filenames));
      } else if(length(ga_lower_bound) != length(ga_stan_filenames)) {
        stop("length of 'ga_lower_bound' should be equal to 1 or the length 
             of 'ga_stan_filenames'");
      }
      names(ga_lower_bound) = names(ga_stan_filenames);
    }
    
    do_these_priors = c(names(hs_stan_filenames), names(ga_stan_filenames));
    
  }
  
  if(include_nonbayes) {
    nonbayes_names = paste0("nonbayes", seq_along(nonbayes_added_weight));
    names(nonbayes_added_weight) = nonbayes_names;
    priors = c(nonbayes_names, do_these_priors);
  } else {
    nonbayes_added_weight = NA;
    priors = do_these_priors;
  }
  
  tiny_positive = .Machine$double.eps;
  priors = factor(priors) %>% fct_inorder
  
  summarized_performance = 
    crossing(priors = priors, 
             crossing(sim_id = 1:n_sim,
                      array_id = array_id,
                      scenario_id = scenario_id, 
                      true_prob_curve_id = true_prob_curve_id, 
                      predictor_dist_id = predictor_dist_id,
                      n_training = n_training,
                      # Were sparse categories observed?
                      obs_sparse_group = NA_real_,
                      # What is probability of sparse categories?
                      prob_sparse_group = NA_real_,
                      #root mean squared error of pointwise probabilities:
                      rmse = NA_real_,
                      #bias using pointwise probabilities:
                      bias = NA_real_,
                      #oos pointwise log-likelihood using model:
                      loglik = NA_real_,
                      #oos pointwise kl divergence:
                      kl_div = NA_real_,
                      #oos log-likelihood using empiric mean of y:
                      loglik_intercept = NA_real_,
                      #oos 50% coverage 
                      coverage_50 = NA_real_));
  
  summarized_bayesian_performance = 
    crossing(priors = priors, 
             crossing(sim_id = 1:n_sim,
                      array_id = array_id,
                      scenario_id = scenario_id, 
                      true_prob_curve_id = true_prob_curve_id, 
                      predictor_dist_id = predictor_dist_id,
                      n_training = n_training,
                      # Were sparse categories observed?
                      obs_sparse_group = NA_real_,
                      # What is probability of sparse categories?
                      prob_sparse_group = NA_real_,
                      #Number of breaks in the predictor used:
                      n_breaks = NA_real_,
                      #Either shape parameter (gamma prior) or scale parameter (horseshoe prior):
                      tuning_param_val = NA_real_,
                      #Bayesian RMSE
                      rmse = NA_real_,
                      #Bayesian bias
                      bias = NA_real_,
                      #Bayesian oos log-likehood
                      loglik = NA_real_,
                      #Bayesian KL divergence:
                      kl_div = NA_real_,
                      #Bayesian coverage
                      coverage_50 = NA_real_, 
                      #how many divergent transitions were there?
                      number_divergences = NA_real_,
                      #what was the value of the gelman-rubin diagnostic?
                      rhat = NA_real_, 
                      #where there NaNs? (another symptom of poor mixing)
                      any_nan = NA_real_,
                      #time required to fit each method
                      run_time_secs = NA_real_,
                      #ratio of run time of fastest versus slowest chains
                      chain_relative_diff = NA_real_)) %>%
    filter(!str_detect(priors, "nonbayes")) %>%
    mutate(priors = fct_drop(priors))
  
  
  if(return_summarized_models) {
    
    prespecified_x_validation = 
      predictor_dist(1e5) %>% unique() %>% sort();
    true_prob_prespecified_validation = 
      true_prob_curve(prespecified_x_validation);
    summarized_models = 
      crossing(priors = priors, 
               crossing(sim_id = 1:n_sim,
                        x = prespecified_x_validation, 
                        array_id = array_id, 
                        scenario_id = scenario_id, 
                        true_prob_curve_id = true_prob_curve_id, 
                        predictor_dist_id = predictor_dist_id,
                        n_training = n_training,
                        true_prob = NA_real_, 
                        fitted_prob = NA_real_,
                        lower50 = NA_real_,
                        upper50 = NA_real_));
  } else {
    summarized_models = NULL;
  }
  
  stan_compiled = F;
  begin_sim = Sys.time();
  i=1;
  
  if(fit_methods) {
    cat("\n######################################\n");
    cat("# the true probability curve is defined by the function: \n");
    cat("# ")
    print(true_prob_curve);
    cat("# the true distribution of the predictor is defined by the function: \n");
    cat("# ")
    print(predictor_dist);
    cat("# the number of observations used for training is", n_training, "\n");
    cat("# the value of 'random_seed' is", random_seed,"\n");
    cat("# the value of 'n_sim' is", n_sim,"\n");
    cat("######################################\n\n");
  }
  
  for(i in 1:n_sim) {
    
    if(fit_methods) {
      
      # Draw data ----
      set.seed(data_seeds[i]);
      x = predictor_dist(n_training);
      length_unique_x = length(unique(x));
      true_prob = true_prob_curve(x);
      if(any(true_prob > 1 | true_prob < 0)) {
        stop(paste0("The function 'true_prob_curve' generated invalid probabilities, i.e. outside [0,1], given the following inputs: ", paste0(formatC(x[which(true_prob > 1 | true_prob < 0)],format = "f", digits = 4), collapse = ", ")));
      }
      y = rbinom(n_training, 1, prob = true_prob);
      mean_y = pmax(pmin(mean(y), 1 - tiny_positive), tiny_positive);
      
      x_validation = predictor_dist(n_validation);
      true_prob_validation = true_prob_curve(x_validation);
      if(any(true_prob_validation > 1 | true_prob_validation < 0)) {
        stop(paste0("The function 'true_prob_curve' generated invalid probabilities, i.e. outside [0,1], given the following inputs: ", paste0(formatC(true_prob_validation[which(true_prob_validation > 1 | true_prob_validation < 0)],format = "f", digits = 4), collapse = ", ")));
      }
      true_prob_validation = pmax(pmin(true_prob_validation, 1 - tiny_positive), tiny_positive)
      
      curr_loglik_intercept = 
        mean(true_prob_validation * log(mean_y) + 
               (1 - true_prob_validation)  * log(1 - mean_y) -
               tiny_positive * log(tiny_positive));
      
      # Non-Bayesian Isotonic Regression first ----
      if(include_nonbayes) {
        for(curr_prior in nonbayes_names) {
          
          cat("\n######################################\n");
          cat("# iteration:",i, "/", n_sim,"\n");
          cat("# the current data seed is",data_seeds[i],"\n");
          cat("#", curr_prior,"\n");
          cat("######################################\n\n");
          
          curr_row_performance = 
            with(summarized_performance, 
                 which(sim_id == i & 
                         priors == curr_prior));
          
          if(return_summarized_models) {
            curr_rows_model = 
              with(summarized_models, 
                   which(sim_id == i & 
                           priors == curr_prior));
          }
          
          
          data_grouped_iso = 
            tibble(x = x,
                   y = y) %>%
            group_by(x) %>%
            summarize(wt = length(y),
                      unweighted_y = mean(y), 
                      weighted_y = (sum(y) + (0.5 * nonbayes_added_weight[curr_prior])) /
                        (length(y) + nonbayes_added_weight[curr_prior])) %>%
            arrange(x)
          
          if(nrow(data_grouped_iso) > 1) {
            # If just one category, reduces to Wilson's confidence interval, i.e.
            # inverted score
            curr_fit <- 
              quickIsotone(doseResponse(y = pull(data_grouped_iso, unweighted_y), 
                                        x = pull(data_grouped_iso, x), 
                                        wt = pull(data_grouped_iso, wt)), 
                           estfun = oldPAVA, 
                           conf = 0.5);
          } else {
            
            foo = binom.confint(x = pull(data_grouped_iso, unweighted_y) * 
                                  pull(data_grouped_iso, wt),
                                n = pull(data_grouped_iso, wt), 
                                conf.level = 0.5, 
                                method = "wilson");
            curr_fit <-
              tibble(x = pull(data_grouped_iso, x), 
                     y = pull(data_grouped_iso, unweighted_y),
                     lower50conf = foo$lower[1],
                     upper50conf = foo$upper[1])
          }
          
          model_performance = 
            validate_isotonic(x_validation = x_validation, 
                              true_prob_validation = true_prob_validation, 
                              fitted_isotonic = curr_fit);
          
          summarized_performance[curr_row_performance, c("rmse","bias","loglik","kl_div","coverage_50")] = 
            model_performance$summary_fitted_isotonic[,c("rmse","bias","loglik","kl_div","coverage")];
          summarized_performance[curr_row_performance, "loglik_intercept"] = 
            curr_loglik_intercept;
          
          #Model summary ----
          if(return_summarized_models) {
            model_summary = 
              validate_isotonic(x = prespecified_x_validation, 
                                true_prob = true_prob_prespecified_validation, 
                                fitted_isotonic = curr_fit);
            
            summarized_models[curr_rows_model,c("x","true_prob","fitted_prob","lower50","upper50")] = 
              model_summary$fitted_isotonic[,c("x_validation","true_prob_validation","y","lower50conf","upper50conf")] %>%
              arrange(x_validation);
            rm(curr_rows_model, model_summary);
          }
          
          rm(data_grouped_iso, model_performance, curr_fit, curr_row_performance, curr_prior);
        }
      }
      if(length(do_these_priors)) {
        
        if(!stan_compiled) {
          # Compile templates ----
          begin_compile = Sys.time();
          
          if(length(hs_stan_filenames)) {
            hs_stan_templates = vector("list",length = length(hs_stan_filenames));
            names(hs_stan_templates) = names(hs_stan_filenames);
            for(j in names(hs_stan_filenames)) {
              hs_stan_templates[[j]] = 
                bayesian_isotonic(data_grouped = tibble(n = 1 + numeric(5), 
                                                        y = numeric(5)),
                                  stan_fit = NA, 
                                  stan_path = hs_stan_filenames[j], 
                                  stan_args = 
                                    list(
                                      local_dof_stan = 1, 
                                      global_dof_stan = 1,
                                      alpha_scale_stan = 1),
                                  sample_from_prior_only = T,
                                  n_mc_warmup = 1e2, 
                                  n_mc_samps = 1e2, 
                                  return_as_stan_object = T,
                                  stan_seed = stan_seeds[i]);
            }
          }
          
          if(length(ga_stan_filenames)) {
            ga_stan_templates = vector("list",length = length(ga_stan_filenames));
            names(ga_stan_templates) = names(ga_stan_filenames);
            for(j in names(ga_stan_filenames)) {
              ga_stan_templates[[j]] =
                bayesian_isotonic(data_grouped = tibble(n = 1 + numeric(5), 
                                                        y = numeric(5)),
                                  stan_fit = NA, 
                                  stan_path = ga_stan_filenames[j], 
                                  stan_args = list(
                                    alpha_shape_stan = 1,
                                    tiny_positive_stan = 1e-5),
                                  sample_from_prior_only = T,
                                  n_mc_warmup = 1e2, 
                                  n_mc_samps = 1e2, 
                                  return_as_stan_object = T,
                                  stan_seed = stan_seeds[i]);
            }
          }
          end_compile = Sys.time();  
          
          
          stan_compiled = T;
        } 
        
        # PHIL: DELETE THESE LINES BELOW ----
        # curr_prior = (do_these_priors)[1];
        # PHIL: DELETE THESE LINES ABOVE ----
        
        data_grouped = 
          tibble(x = x, 
                 y = y) %>%
          group_by(x) %>%
          summarize(n = length(x),
                    y = sum(y)) %>%
          arrange(x) %>%
          mutate(x_cat = 1:n())
        
        validation_data_ungrouped <- 
          left_join(
            tibble(x = x_validation) %>%
              mutate(row_num = 1:n()),
            data_grouped %>%
              select(x, x_cat)) %>%
          arrange(x) %>%
          fill(x_cat, 
               .direction = "downup") %>%
          arrange(row_num)
        
        prespecified_validation_data_ungrouped <- 
          left_join(
            tibble(x = prespecified_x_validation) %>%
              mutate(row_num = 1:n()),
            data_grouped %>%
              select(x, x_cat), 
            by = "x") %>%
          arrange(x) %>%
          fill(x_cat, 
               .direction = "downup") %>%
          arrange(row_num) 
        
        for(curr_prior in do_these_priors) {
          
          cat("\n######################################\n");
          cat("# iteration:",i, "/", n_sim,"\n");
          cat("# the current data seed is",data_seeds[i],"\n");
          cat("# the current stan seed is",stan_seeds[i],"\n");
          cat("#", curr_prior, "prior","\n");
          cat("######################################\n\n");
          
          curr_row_performance = 
            with(summarized_performance, 
                 which(sim_id == i & 
                         priors == curr_prior));
          
          curr_row_bayesian_performance = 
            with(summarized_bayesian_performance, 
                 which(sim_id == i & 
                         priors == curr_prior));
          
          if(return_summarized_models) {
            curr_rows_model = 
              with(summarized_models, 
                   which(sim_id == i & 
                           priors == curr_prior));
          }
          
          
          if(curr_prior %in% names(hs_stan_filenames)) {
            # horseshoe
            stan_args = 
              list(
                local_dof_stan = 1, 
                global_dof_stan = 1,
                alpha_scale_stan = 
                  solve_for_hs_scale(
                    target_mean = 0.5 / length_unique_x,
                    local_dof = 1, 
                    global_dof = 1, 
                    slab_precision = 1,
                    n = n_training,
                    sigma = 2
                  )$scale
              );
            
            summarized_bayesian_performance[curr_row_bayesian_performance, "tuning_param_val"] = 
              stan_args$alpha_scale_stan;
            
            stan_template = hs_stan_templates[[curr_prior]];
          } else {
            # gamma
            stan_args = 
              list(
                alpha_shape_stan = 0.5 / length_unique_x,
                tiny_positive_stan = ga_lower_bound[curr_prior]);
            
            summarized_bayesian_performance[curr_row_bayesian_performance, "tuning_param_val"] = 
              stan_args$alpha_shape_stan;
            
            stan_template = ga_stan_templates[[curr_prior]];
            
          }
          
          summarized_performance[curr_row_performance, "obs_sparse_group"] = 
            summarized_bayesian_performance[curr_row_bayesian_performance, "obs_sparse_group"] =
            any(pull(data_grouped, y) == 0 | 
                  pull(data_grouped, y) == pull(data_grouped, n))
          
          summarized_performance[curr_row_performance, "prob_sparse_group"] =  
            summarized_bayesian_performance[curr_row_bayesian_performance, "prob_sparse_group"] =
            tibble(x = x, 
                   true_prob = true_prob) %>%
            group_by(x) %>%
            summarize(prob_homogenous = 
                        prod(true_prob) + prod(1 - true_prob)) %>%
            summarize(prob_any_homogenous = 1 - prod(1 - prob_homogenous)) %>%
            pull(prob_any_homogenous)
          
          #Fit Stan model ----
          curr_fit = bayesian_isotonic(data_grouped = data_grouped,
                                       stan_fit = stan_template,
                                       stan_path = NA, 
                                       stan_args = stan_args, 
                                       # 'conf_level' must remain hard coded for this simulator
                                       # because the function 'validate_bayesian_isotonic' assumes 
                                       # this level
                                       conf_level = 0.50, 
                                       conf_level_direction = "both",
                                       sample_from_prior_only = F,
                                       n_mc_warmup = n_mc_warmup, 
                                       n_mc_samps = n_mc_samps, 
                                       mc_chains = mc_chains, 
                                       mc_thin = mc_thin, 
                                       mc_stepsize = mc_stepsize, 
                                       mc_adapt_delta = mc_adapt_delta,
                                       mc_max_treedepth = mc_max_treedepth,
                                       verbose = T,
                                       stan_seed = stan_seeds[i]);
          
          #Non-Bayesian Model performance ----
          
          fitted_isotonic =
            bind_cols(x = pull(data_grouped, x), 
                      y = apply(curr_fit$all_draws$xi, 2, quantile, p = 0.50, na.rm = TRUE),  
                      lower50conf = apply(curr_fit$all_draws$xi, 2, quantile, p = 0.25, na.rm = TRUE), 
                      upper50conf = apply(curr_fit$all_draws$xi, 2, quantile, p = 0.75, na.rm = TRUE))              
          
          model_performance = 
            validate_isotonic(x_validation = validation_data_ungrouped$x, 
                              true_prob_validation = true_prob_validation, 
                              fitted_isotonic = fitted_isotonic);
          
          summarized_performance[curr_row_performance, 
                                 c("rmse","bias","loglik","kl_div","coverage_50")] = 
            model_performance$summary_fitted_isotonic[,c("rmse","bias","loglik","kl_div","coverage")];
          summarized_performance[curr_row_performance, "loglik_intercept"] = 
            curr_loglik_intercept;
          rm(fitted_isotonic, curr_row_performance, model_performance);
          
          # Bayesian model performance
          bayesian_model_performance = 
            validate_bayesian_isotonic(x_validation_category = validation_data_ungrouped$x_cat,
                                       true_prob_validation = true_prob_validation, 
                                       draws_xi = curr_fit$all_draws$xi);
          
          
          summarized_bayesian_performance[curr_row_bayesian_performance, 
                                          c("rmse","bias","loglik","kl_div","coverage_50")] = 
            bayesian_model_performance$summary_fitted_isotonic[c("rmse","bias","loglik","kl_div","coverage")];
          
          
          summarized_bayesian_performance[curr_row_bayesian_performance, "number_divergences"] =
            curr_fit$number_divergences;
          summarized_bayesian_performance[curr_row_bayesian_performance, "rhat"] = curr_fit$max_rhat;
          summarized_bayesian_performance[curr_row_bayesian_performance, "any_nan"] =
            curr_fit$any_nan;
          summarized_bayesian_performance[curr_row_bayesian_performance, "run_time_secs"] = 
            curr_fit$total_run_time_secs;
          summarized_bayesian_performance[curr_row_bayesian_performance, "chain_relative_diff"] = 
            min(curr_fit$chain_run_times_secs) / max(curr_fit$chain_run_times_secs);
          
          #Model summary ----
          if(return_summarized_models) {
            model_summary = 
              validate_bayesian_isotonic(
                x = prespecified_validation_data_ungrouped$x_cat, 
                true_prob = true_prob_prespecified_validation, 
                draws_xi = curr_fit$all_draws$xi);
            
            summarized_models[curr_rows_model,
                              c("x","true_prob","fitted_prob","lower50","upper50")] = 
              bind_cols(model_summary$returned_predictions[,c("x_cat","true_prob","fitted_prob","lower50","upper50")],
                        x = prespecified_validation_data_ungrouped$x) %>%
              select(x, true_prob, fitted_prob, lower50, upper50) %>%
              arrange(x);
            rm(curr_rows_model, model_summary);
          }
          rm(bayesian_model_performance, curr_row_bayesian_performance);
          
          rm(curr_fit, stan_args, stan_template);
        }
        rm(curr_prior);
      }
      rm(data_grouped, x, true_prob, y, mean_y, 
         prespecified_validation_data_ungrouped,
         x_validation, validation_data_ungrouped, 
         true_prob_validation, curr_loglik_intercept);
      
    } else if(i == 1) {
      cat("Conducting dry run only. No simulations were run.\n");
    }
  }
  
  # +++ Report out results ----  
  
  list(summarized_performance = summarized_performance,
       summarized_bayesian_performance = summarized_bayesian_performance,
       summarized_models = summarized_models, 
       nonbayes_added_weight = nonbayes_added_weight,
       hs_stan_filenames = hs_stan_filenames, 
       ga_stan_filenames = hs_stan_filenames,
       ga_lower_bound = ga_lower_bound,
       n_training = n_training,
       n_validation = n_validation,
       true_prob_curve = true_prob_curve,
       predictor_dist = predictor_dist, 
       random_seed = random_seed,
       data_seeds = data_seeds,
       stan_seeds = stan_seeds,
       total_run_time_secs = difftime(Sys.time(), begin_all, units = "secs"));
} 


