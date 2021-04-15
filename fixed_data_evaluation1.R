# code for Section 3.1: Fixed data evaluation
# set run_sims = T to conduct Fixed-data evaluation 1 in Section 3.1
# set run_sims = F if you have already run the numerical study and want to 
# process the results (the files must be in the same place that they were
# written to, namely the local folder called 'out'). The last call in the 
# script generates a latex table that is Table 1 in the manuscript
run_sims = T;

library(rstan);
library(tidyverse);
library(cir);
library(glue);
library(kableExtra);

#Helpful to print warnings when they occur for debugging
options(warn = 1);

#Recommended options from rstan:
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

source("1functions_methods.R");

data_features = 
  tibble(n_cats = 1,
         n_per_cat = c(50, 500))


n_mc_warmup = 2.5e3;#
n_mc_samps = 2.5e3;#

all_data_grouped <- 
  hs_scales <- 
  ga_shapes <- 
  ga_lower_bounds <- 
  vector("list")

for(j in 1:nrow(data_features)) {
  
  curr_n_cats = slice(data_features, j) %>% pull(n_cats)
  curr_n_per_cat = slice(data_features, j) %>% pull(n_per_cat)
  
  hs_scales[[j]] =
    c(horseshoe1 = 1 / (curr_n_cats + 1), 
      horseshoe2 = (1/4) / (curr_n_cats + 1), 
      horseshoe3 = (1/16) / (curr_n_cats + 1), 
      horseshoe4 = (1/64) / (curr_n_cats + 1));
  
  ga_shapes[[j]] =
    c(gamma1 = 1 / (curr_n_cats + 1), 
      gamma2 = (1/4) / (curr_n_cats + 1), 
      gamma3 = (1/16) / (curr_n_cats + 1), 
      gamma4 = (1/64) / (curr_n_cats + 1));
  
  
  ga_lower_bounds[[j]] = 
    c(gamma1 = 0,
      gamma2 = 0,
      gamma3 = 0,
      gamma4 = 0);
  
  
  all_data_grouped[[j]] = 
    tibble(dataset_label = j, 
           n_cats = curr_n_cats, 
           x = 1:curr_n_cats, 
           n = curr_n_per_cat,
           y = round(rep(curr_n_per_cat/2, curr_n_cats))) %>%
    arrange(x) %>%
    mutate(x_cat = 1:n());
  
}
rm(j, curr_n_cats, curr_n_per_cat);

all_data_grouped <- bind_rows(all_data_grouped)


if(run_sims) {
  
  set.seed(1);
  stan_seeds = sample.int(.Machine$integer.max, 1);
  
  hs_stan_filenames = 
    c(horseshoe1 = "iso_horseshoe.stan",
      horseshoe2 = "iso_horseshoe.stan",
      horseshoe3 = "iso_horseshoe.stan",
      horseshoe4 = "iso_horseshoe.stan"); 
  
  ga_stan_filenames = 
    c(gamma1 = "iso_gamma.stan",
      gamma2 = "iso_gamma.stan",
      gamma3 = "iso_gamma.stan",
      gamma4 = "iso_gamma.stan");
  
  do_these_priors = c(names(hs_stan_filenames),
                      names(ga_stan_filenames));
  
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
                              alpha_scale_stan = 1e-2, 
                              slab_precision_stan = 1),
                          sample_from_prior_only = T,
                          n_mc_warmup = 1e2, 
                          n_mc_samps = 1e2, 
                          return_as_stan_object = T);
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
                          return_as_stan_object = T);
    }
  }
  end_compile = Sys.time();  
  
  summarized_performance = 
    crossing(
      dataset_label = pull(all_data_grouped, dataset_label) %>% unique(),
      priors = factor(do_these_priors) %>% fct_inorder(), 
      #Either shape parameter (gamma prior) or scale parameter (horseshoe prior):
      tuning_param_val = NA_real_,
      #What is the left of the support truncated to? (0 means no truncation)
      lower_truncation = NA_real_,
      #how many divergent transitions were there?
      number_divergences = NA_real_,
      #what was the value of the gelman-rubin diagnostic?
      rhat = NA_real_, 
      #how many NaNs were sampled? (another symptom of poor mixing)
      any_nan = NA_real_,
      #time required to fit each method
      run_time_secs = NA_real_,
      #ratio of run time of fastest versus slowest chains
      chain_relative_diff = NA_real_);
  
  alpha_nans = 
    crossing(
      all_data_grouped %>% select(dataset_label, x_cat),
      priors = factor(do_these_priors) %>% fct_inorder(), 
      value = NA_real_) %>%
    pivot_wider(id_cols = c(dataset_label, priors),
                names_from = x_cat, 
                names_prefix = "alpha", 
                values_from = value,
                values_fill = NA_real_);
  
  xi_nans = 
    crossing(
      all_data_grouped %>% select(dataset_label, x_cat),
      priors = factor(do_these_priors) %>% fct_inorder(), 
      value = NA_real_) %>%
    pivot_wider(id_cols = c(dataset_label, priors),
                names_from = x_cat, 
                names_prefix = "xi", 
                values_from = value,
                values_fill = NA_real_);
  
  
  summarized_models = 
    crossing(
      bind_rows(
        all_data_grouped,
        all_data_grouped %>% 
          group_by(dataset_label) %>%
          slice(n()) %>% 
          mutate(x = x + 1, y= NA, x_cat = x_cat + 1)
      ), 
      priors = factor(do_these_priors) %>% fct_inorder(), 
      tuning_param_val = NA_real_,
      lower_truncation = NA_real_,
      #
      alphaj = NA_real_,
      alphaj_lower = NA_real_,
      alphaj_upper = NA_real_,
      xij = NA_real_,
      xij_lower = NA_real_,
      xij_upper = NA_real_) %>%
    arrange(dataset_label, priors, x)
  
  
  for(j in pull(all_data_grouped, dataset_label) %>% unique()) {
    
    curr_n_cats = 
      data_features %>% 
      slice(j) %>%
      pull(n_cats)
    
    curr_n_per_cat = 
      data_features %>% 
      slice(j) %>%
      pull(n_per_cat)
    
    data_grouped = 
      all_data_grouped %>%
      filter(dataset_label == j) %>%
      select(-dataset_label, -n_cats)
    
    curr_hs_scale = hs_scales[[j]]
    curr_ga_shape = ga_shapes[[j]]
    curr_ga_lower_bound = ga_lower_bounds[[j]]
    
    for(curr_prior in do_these_priors) {
      
      cat("\n######################################\n");
      cat("#", curr_prior, "prior :: dataset =", j,  "\n");
      cat("######################################\n\n");
      
      
      if(curr_prior %in% names(hs_stan_filenames)) {
        
        # horseshoe
        stan_args = 
          list(local_dof_stan = 1, 
               global_dof_stan = 1,
               alpha_scale_stan = curr_hs_scale[[curr_prior]],
               slab_precision_stan = 1);
        
        stan_template = hs_stan_templates[[curr_prior]];
      } else {
        
        # gamma
        stan_args = 
          list(
            alpha_shape_stan = curr_ga_shape[[curr_prior]],
            tiny_positive_stan = curr_ga_lower_bound[[curr_prior]]);
        
        stan_template = ga_stan_templates[[curr_prior]];
        
      }
      
      
      curr_row_performance = 
        with(summarized_performance, 
             which(dataset_label == j &
                     priors == curr_prior));
      
      curr_row_nans = 
        with(alpha_nans, 
             which(dataset_label == j &
                     priors == curr_prior));
      
      stopifnot(length(curr_row_performance) == 1)
      
      curr_rows_models <- 
        with(summarized_models,
             which(dataset_label == j &
                     priors == curr_prior));
      
      if(curr_prior %in% names(hs_stan_filenames)) {
        
        summarized_performance[curr_row_performance, "tuning_param_val"] = 
          curr_hs_scale[[curr_prior]];
        summarized_performance[curr_row_performance, "lower_truncation"] = 
          0;
        
        summarized_models[curr_rows_models, "tuning_param_val"] = 
          curr_hs_scale[[curr_prior]];
        summarized_models[curr_rows_models, "lower_truncation"] = 
          0;
      } else {
        summarized_performance[curr_row_performance, "tuning_param_val"] = 
          curr_ga_shape[[curr_prior]];
        summarized_performance[curr_row_performance, "lower_truncation"] = 
          curr_ga_lower_bound[[curr_prior]];
        
        summarized_models[curr_rows_models, "tuning_param_val"] = 
          curr_ga_shape[[curr_prior]];
        summarized_models[curr_rows_models, "lower_truncation"] = 
          curr_ga_lower_bound[[curr_prior]];
      }
      
      # ++ Fit Stan model ----
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
                                   mc_chains = 2, 
                                   mc_thin = 1, 
                                   mc_stepsize = 0.1, 
                                   # +++ phil check below (should be 0.99 for actual run)----
                                   mc_adapt_delta = 0.99,
                                   # +++ phil check above (should be 0.99 for actual run)----
                                   mc_max_treedepth = 15,
                                   verbose = T, 
                                   stan_seed = stan_seeds,
                                   # +++ phil comment below ----
                                   #return_as_stan_object = TRUE
                                   # +++ phil comment above ----
      );
      
      
      # ++ Model performance ----
      
      summarized_performance[curr_row_performance, "number_divergences"] =
        curr_fit$number_divergences;
      summarized_performance[curr_row_performance, "rhat"] = 
        curr_fit$max_rhat;
      summarized_performance[curr_row_performance, "any_nan"] =
        curr_fit$any_nan;
      summarized_performance[curr_row_performance, "run_time_secs"] = 
        curr_fit$total_run_time_secs;
      summarized_performance[curr_row_performance, "chain_relative_diff"] = 
        min(curr_fit$chain_run_times_secs) / max(curr_fit$chain_run_times_secs);
      
      alpha_nans <- 
        alpha_nans %>%
        mutate_at(vars(paste0("alpha",1:(curr_n_cats))),
                  ~ ifelse(dataset_label == j & 
                             priors == curr_prior, 
                           curr_fit$alpha_number_nan[1:(curr_n_cats)], .)) 
      
      xi_nans <- 
        xi_nans %>%
        mutate_at(vars(paste0("xi",1:(curr_n_cats))),
                  ~ ifelse(dataset_label == j & 
                             priors == curr_prior,
                           curr_fit$xi_number_nan[1:(curr_n_cats)], .)) 
      
      summarized_models[curr_rows_models, "xij"] = 
        c(apply(curr_fit$all_draws$xi, 2, median), 1)
      summarized_models[curr_rows_models, "xij_lower"] = 
        c(apply(curr_fit$all_draws$xi, 2, quantile, p = 0.025), 1)
      summarized_models[curr_rows_models, "xij_upper"] = 
        c(apply(curr_fit$all_draws$xi, 2, quantile, p = 0.975), 1)
      
      summarized_models[curr_rows_models, "alphaj"] = 
        apply(curr_fit$all_draws$alpha, 2, median)
      summarized_models[curr_rows_models, "alphaj_lower"] = 
        apply(curr_fit$all_draws$alpha, 2, quantile, p = 0.025)
      summarized_models[curr_rows_models, "alphaj_upper"] = 
        apply(curr_fit$all_draws$alpha, 2, quantile, p = 0.975)
      
      
    } 
    
    write_csv(filter(summarized_performance, 
                     dataset_label <= j),
              path = "out/exemplar1_fixed_data_performance.csv",
              append = FALSE);
    
    write_csv(filter(summarized_models,
                     dataset_label <= j),
              path = "out/exemplar1_fixed_data_models.csv",
              append = FALSE);
    
    write_csv(filter(xi_nans, 
                     dataset_label <= j),
              path = "out/exemplar1_xi_nans.csv",
              append = FALSE);
    
    write_csv(filter(alpha_nans, 
                     dataset_label <= j),
              path = "out/exemplar1_alpha_nans.csv",
              append = FALSE);
    
    rm(stan_template, stan_args, curr_prior);
    rm(curr_row_performance, curr_row_nans, curr_rows_models);
    
    
  }
  rm(curr_ga_shape, curr_hs_scale, curr_ga_lower_bound);
  rm(data_grouped, curr_n_cats, curr_n_per_cat)
  
} else {
  
  summarized_performance = 
    read_csv("out/exemplar1_fixed_data_performance.csv")
  
  summarized_models = 
    read_csv("out/exemplar1_fixed_data_models.csv")
  
  xi_nans = 
    read_csv("out/exemplar1_xi_nans.csv")
  
  alpha_nans =
    read_csv("out/exemplar1_alpha_nans.csv")
  
  
  summarized_performance %>% 
    group_by(dataset_label, priors) %>% 
    summarize(median_div = median(number_divergences),
              q1_div = quantile(number_divergences, p = 0.25),
              q3_div = quantile(number_divergences, p = 0.75),
              max_div = max(number_divergences),
              mean_any_nan = mean(any_nan > 0),
              median_rhat = median(rhat, na.rm = T),
              q1_rhat = quantile(rhat, p = 0.25),
              q3_rhat = quantile(rhat, p = 0.75),
              max_rhat = max(rhat),
              median_runtime = median(run_time_secs),
              q1_runtime = quantile(run_time_secs, p = 0.25),
              q3_runtime = quantile(run_time_secs, p = 0.75),
              max_runtime = max(run_time_secs)) %>%
    ungroup() %>% 
    mutate(
      max_rhat = ifelse(max_rhat < 1e9, formatC(max_rhat, format = "f", digits = 2), "$>10^9$"),
      div = glue("{round(median_div)}({max_div})"), 
      nan = formatC(mean_any_nan, format = "f", digits = 2), 
      rhat = glue("{formatC(median_rhat,format = 'f', digits = 2)}({max_rhat})"),
      runtime = glue("{round(median_runtime)}({round(max_runtime)})")) %>%
    select(dataset_label, priors, div:runtime)
  
  
  
  summarized_models <-
    summarized_models %>%
    mutate(
      dataset_pretty_label = 
        glue("list(K == {n_cats}, n[j] == {n})") %>% 
        as.character() %>%
        factor() %>%
        fct_inorder(),
      priors_pretty_label = 
        case_when(
          priors == "gamma1" ~ "GA[1]",
          priors == "gamma2" ~ "GA[2]",
          priors == "gamma3" ~ "GA[3]",
          priors == "gamma4" ~ "GA[4]",
          priors == "horseshoe1" ~ "HS",
        ),
      tuning_param_val_pretty = 
        case_when(
          str_detect(priors, "gamma") ~ 
            glue("$s = 1/{1 / tuning_param_val}$"),
          str_detect(priors, "horseshoe") ~ 
            glue("$c = 1/{1 / tuning_param_val}$"),
        ),
      alphaj_pretty = 
        formatC(alphaj, format = 'g', digits = 2),
      #glue("{formatC(alphaj, format = 'g', digits = 2)}({formatC(alphaj_lower, format = 'g', digits = 2)},{formatC(alphaj_upper, format = 'g', digits = 2)})"),
      xij_pretty = 
        glue("{formatC(xij, format = 'f', digits = 2)}({formatC(xij_lower, format = 'f', digits = 2)},{formatC(xij_upper, format = 'f', digits = 2)})"))
  
  
  
  full_join(
    summarized_models %>% 
      select(n, priors, tuning_param_val, tuning_param_val_pretty, x_cat, alphaj_pretty) %>%
      pivot_wider(names_from = x_cat, 
                  names_prefix = "alpha",
                  values_from = alphaj_pretty),
    summarized_models %>% 
      filter(x_cat == 1) %>%
      select(n, priors, tuning_param_val_pretty, xij_pretty),
    by = c("n", "priors", "tuning_param_val_pretty")) %>% 
    mutate(priors = 
             case_when(
               str_detect(priors, "horseshoe") ~ "HSIPV",
               str_detect(priors, "gamma") ~ "GAIPV",
             ) %>%
             factor(levels = c("GAIPV", "HSIPV"))) %>%
    arrange(n, priors, desc(tuning_param_val)) %>%
    select(-tuning_param_val) %>%
    rename(`$n$` = n, 
           Prior = priors, 
           `Tuning parameter` = tuning_param_val_pretty,
           `$\\alpha_1$` = alpha1,
           `$\\alpha_2$` = alpha2,
           `$\\xi_1$ (95\\% CI)` = xij_pretty) %>%
    knitr::kable(format = "latex",
                 booktabs = T,
                 longtable = F,
                 escape = F,
                 linesep = c("", "", "", "\\addlinespace")) %>%
    kable_styling(latex_options = c("HOLD_position"),
                  full_width = F,
                  font_size = 11) %>%
    add_header_above(c(" " = 3, "Posterior median" = 3))
}

