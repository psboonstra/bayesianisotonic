
make_grouped_data <- function(x, y, breaks = NULL) {
  
  if(is.null(breaks)) {
    message("Categorizing 'x' into five groups. Consider choosing your own value.")
    breaks <- 
      quantile(x, probs = seq(0, 1, length = 6)) %>%
      as.numeric() %>%
      replace(which.min(.), -Inf) %>% 
      replace(which.max(.), Inf)
  }
  
  tibble(x = cut(x, breaks = breaks, right = F), 
         y = y) %>%
    group_by(x) %>%
    summarize(n = length(x),
              y = sum(y)) %>%
    arrange(x) %>%
    mutate(x_cat = 1:n())   
}

bayesian_isotonic = function(data_grouped = NULL,
                             stan_fit = NA, 
                             stan_path = NA,
                             stan_args = list(
                               local_dof_stan = 1, 
                               global_dof_stan = 1,
                               alpha_scale_stan = 1),
                             sample_from_prior_only = F,
                             conf_level = 0.50, 
                             conf_level_direction = "both",
                             sig_threshold = c(0.005, 0.01, 0.05),
                             verbose = F,
                             n_mc_warmup = 2.5e3, 
                             n_mc_samps = 5e3, 
                             mc_chains = 1, 
                             mc_thin = 1, 
                             mc_stepsize = 0.1, 
                             mc_adapt_delta = 0.99,
                             mc_max_treedepth = 15,
                             return_as_stan_object = F,
                             tol = .Machine$double.eps^0.5,
                             stan_seed = sample.int(.Machine$integer.max, 1)) {
  
  require(tidyverse);require(rstan);
  
  stopifnot("data.frame" %in% class(data_grouped))
  stopifnot(c("y","n") %in% colnames(data_grouped));
  stopifnot(all(pull(data_grouped,y) >= 0) && 
              all(pull(data_grouped,y) <= pull(data_grouped,n)));
  
  curr_fit <- stan(file = stan_path,
                   fit = stan_fit,
                   data = c(list(n_groups_stan = nrow(data_grouped),
                                 n_per_group_stan = as.array(pull(data_grouped,n)),
                                 y_stan = as.array(pull(data_grouped,y)),
                                 only_prior_stan = as.integer(sample_from_prior_only)),
                            stan_args), 
                   warmup = n_mc_warmup, 
                   iter = n_mc_samps + n_mc_warmup, 
                   chains = mc_chains, 
                   thin = mc_thin,
                   seed = stan_seed,
                   verbose = F,
                   control = list(stepsize = mc_stepsize,
                                  adapt_delta = mc_adapt_delta,
                                  max_treedepth = mc_max_treedepth));
  
  number_divergences = count_stan_divergences(curr_fit);
  max_rhat = max(summary(curr_fit)$summary[,"Rhat"], na.rm = TRUE)
  
  if(!return_as_stan_object) {
    chain_run_times_secs = rowSums(rstan::get_elapsed_time(curr_fit));
    total_run_time_secs = max(chain_run_times_secs);
    foo = rstan::extract(curr_fit);
    
    xi_number_nan = colSums(is.na(foo$xi))
    alpha_number_nan = colSums(is.na(foo$alpha))
    mean_prob = colMeans(foo$xi, na.rm = T);
    
    if(conf_level_direction == "both") {
      quantile_probs <- 
        apply(foo$xi, 
              2, 
              quantile,
              p = 1/2 + c(-conf_level, 0, conf_level)/2, na.rm = T);
    } else if(conf_level_direction == "lower") {
      quantile_probs <- 
        rbind(apply(foo$xi,
                    2,
                    quantile,
                    p = c(1 - conf_level, 1/2), na.rm = T),
              "100%" = 1);
    } else {
      quantile_probs <- 
        rbind("0%" = 0, 
              apply(foo$xi,
                    2,
                    quantile,
                    p = c(1/2, conf_level), na.rm = T));
    }
    quantile_probs = t(quantile_probs);
    colnames(quantile_probs) = c("model_lower_ci_prob", "model_median_prob", "model_upper_ci_prob");
  } else {
    foo = rstan::extract(curr_fit);
    xi_number_nan = colSums(is.na(foo$xi))
    alpha_number_nan = colSums(is.na(foo$alpha))
  }
  
  
  
  if(number_divergences > 0) {
    warning(paste0("there were ", number_divergences, " divergent transitions"));
  }
  if(any(xi_number_nan > 0) | any(alpha_number_nan > 0))  {
    warning(paste0("there were ", max(c(xi_number_nan,alpha_number_nan)), " draws in which one or more elements of xi were NaN"));
  }
  
  if(return_as_stan_object) {
    curr_fit;
  } else {
    
    data_grouped = 
      data_grouped %>%
      mutate(emp_mean_prob = y/n) %>%
      bind_cols(model_mean_prob = mean_prob,
                as_tibble(quantile_probs));
    
    draws_delta = t(apply(cbind(0,0,foo$xi),1,diff))[,-1,drop = F];
    for(i in 1:length(sig_threshold)) {
      data_grouped = 
        bind_cols(data_grouped,
                  !!sym(paste0("prob_delta_gt_",sig_threshold[i])) := colMeans(draws_delta > sig_threshold[i], na.rm = T));
    }
    
    if(verbose) {
      c(list(data_grouped = data_grouped,
             conf_level = conf_level),
        stan_args,
        list(sample_from_prior_only = sample_from_prior_only,
             number_divergences = number_divergences,
             max_rhat = max_rhat,
             xi_number_nan = xi_number_nan,
             alpha_number_nan = alpha_number_nan,
             any_nan = max(c(xi_number_nan, alpha_number_nan) > 0),
             all_draws = foo,
             chain_run_times_secs = chain_run_times_secs,
             total_run_time_secs = total_run_time_secs));
    } else {
      c(list(data_grouped = data_grouped,
             conf_level = conf_level),
        stan_args,
        list(sample_from_prior_only = sample_from_prior_only,
             number_divergences = number_divergences,
             max_rhat = max_rhat,
             xi_number_nan = xi_number_nan,
             alpha_number_nan = alpha_number_nan,
             any_nan = max(c(xi_number_nan, alpha_number_nan) > 0),
             all_draws = NA,
             chain_run_times_secs = chain_run_times_secs,
             total_run_time_secs = total_run_time_secs));
      
    }
  }
}



count_stan_divergences = function(stan_fit) {
  foo = get_sampler_params(stan_fit, inc_warmup = FALSE);
  n_draws = lapply(foo, nrow)[[1]];
  sum(unlist(lapply(foo,"[",i = 1:n_draws, j = "divergent__")));
}


solve_for_hs_scale = function(target_mean,
                               local_dof = 1,
                               global_dof = 1,
                               slab_precision = 1,
                               n,
                               sigma = 2,
                               tol = .Machine$double.eps,
                               max_iter = 100, 
                               n_sim = 1e6
) {
  
  do_local = (local_dof > 0);
  do_global = (global_dof > 0);
  if(do_local) {
    lambda = rt(n_sim,df = local_dof);
  } else {
    lambda = rep(1, n_sim);
  }
  if(do_global) {
    lambda = lambda * rt(n_sim,df = global_dof);
  } 
  abs_normal_draws = abs(rnorm(n_sim));
  stopifnot(target_mean > 0 && target_mean < 1);#Ensure proper bounds
  log_scale = diff_target = numeric(max_iter);
  log_scale[1] = log(target_mean/(1 - target_mean)*sigma/sqrt(n));
  prior_scales = 1 / sqrt(slab_precision + 1/(exp(2*log_scale[1]) * lambda^2));
  kappa = 1/(1+n*prior_scales/sigma^2);
  diff_target[1] = mean(1-kappa) - target_mean;
  log_scale[2] = 0.5 + log_scale[1];
  prior_scales = 1 / sqrt(slab_precision + 1/(exp(2*log_scale[2]) * lambda^2));
  kappa = 1/(1+n*prior_scales/sigma^2);
  diff_target[2] = mean(1-kappa) - target_mean;
  i=2;
  while(T) {
    i = i+1;
    if(i > max_iter) {i = i-1; break;}
    log_scale[i] =
      log_scale[i-1] - 
      diff_target[i-1]*(log_scale[i-1]-log_scale[i-2])/(diff_target[i-1]-diff_target[i-2]);
    prior_scales =
      1 / sqrt(slab_precision + 1/(exp(2*log_scale[i]) * lambda^2));
    kappa = 1/(1+n*prior_scales/sigma^2);
    diff_target[i] = mean(1-kappa) - target_mean;
    if(abs(diff_target[i]-diff_target[i-1]) < tol) {break;}
  }
  
  list(scale = exp(log_scale[i]),
       achieved_mean = mean(1-kappa),
       target_mean = target_mean,
       diff_from_target = abs(diff_target[i]),
       iter = i);
}

