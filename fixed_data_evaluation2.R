# code for Section 3.2: Fixed data evaluation 2
# set run_sims = T to conduct Fixed-data evaluation 2 in Section 3.2
# set run_sims = F if you have already run the numerical study and want to 
# process the results (the files must be in the same place that they were
# written to, namely the local folder called 'out'). The last call in the 
# script generates Figures 1 and 2 in the manuscript
run_sims = T;

library(rstan);
library(tidyverse);
library(cir);
library(glue);

#Helpful to print warnings when they occur for debugging
options(warn = 1);

#Recommended options from rstan:
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

source("functions.R");

data_features = 
  tibble(n_cats = c(10, 5, 10), 
         n_per_cat = c(8, 16, 32))

n_sim = 50;#
n_mc_warmup = 2.5e3;#
n_mc_samps = 2.5e3;#

all_data_grouped <- 
  hs_scales <- 
  ga_shapes <- 
  dir_shapes <- 
  ga_lower_bounds <- 
  vector("list")

for(j in 1:nrow(data_features)) {
  
  curr_n_cats = slice(data_features, j) %>% pull(n_cats)
  curr_n_per_cat = slice(data_features, j) %>% pull(n_per_cat)
  
  hs_scales[[j]] =
    c(horseshoe1 =
        solve_for_hs_scale1(
          target_mean = 0.5 / (curr_n_cats + 1),
          local_dof = 1, 
          global_dof = 1, 
          slab_precision = 1,
          n = curr_n_per_cat * curr_n_cats,
          sigma = 2
        )$scale);
  
  dir_shapes[[j]] = 
    c(dirichlet1 = 0.5 / (curr_n_cats + 1));
  
  ga_shapes[[j]] =
    c(gamma1 = 0.5 / (curr_n_cats + 1), 
      gamma2 = 0.5 / (curr_n_cats + 1), 
      gamma3 = 0.5 / (curr_n_cats + 1), 
      gamma4 = 0.5 / (curr_n_cats + 1));
  
  
  ga_lower_bounds[[j]] = 
    c(gamma1 = .Machine$double.eps,
      gamma2 = 10^(log10(.Machine$double.eps) - 1), 
      gamma3 = 10^(log10(.Machine$double.eps) - 2), 
      gamma4 = 0);
  
  
  all_data_grouped[[j]] = 
    tibble(dataset_label = j, 
           n_cats = curr_n_cats, 
           x = 1:curr_n_cats, 
           n = curr_n_per_cat,
           y = round(seq(0, curr_n_per_cat, length = curr_n_cats))) %>%
    arrange(x) %>%
    mutate(x_cat = 1:n());
  
}
rm(j, curr_n_cats, curr_n_per_cat);

all_data_grouped <- bind_rows(all_data_grouped)

if(run_sims) {
  
  set.seed(1);
  stan_seeds = sample.int(.Machine$integer.max, n_sim);
  
  hs_stan_filenames = 
    c(horseshoe1 = "iso_horseshoe.stan"); 
  
  ga_stan_filenames = 
    c(gamma1 = "iso_gamma.stan",
      gamma2 = "iso_gamma.stan",
      gamma3 = "iso_gamma.stan",
      gamma4 = "iso_gamma.stan");
  
  
  dir_stan_filenames = 
    c(dirichlet1 = "iso_dirichlet.stan");
  
  do_these_priors = c(names(hs_stan_filenames),
                      names(ga_stan_filenames),
                      names(dir_stan_filenames)
  );
  
  
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
  
  if(length(dir_stan_filenames)) {
    dir_stan_templates = vector("list",length = length(dir_stan_filenames));
    names(dir_stan_templates) = names(dir_stan_filenames);
    for(j in names(dir_stan_filenames)) {
      dir_stan_templates[[j]] =
        bayesian_isotonic(data_grouped = tibble(n = 1 + numeric(5), 
                                                y = numeric(5)),
                          stan_fit = NA, 
                          stan_path = dir_stan_filenames[j], 
                          stan_args = list(
                            alpha_shape_stan = 1),
                          sample_from_prior_only = T,
                          n_mc_warmup = 1e2, 
                          n_mc_samps = 1e2, 
                          return_as_stan_object = T);
    }
  }
  end_compile = Sys.time();  
  
  
  
  summarized_performance = 
    crossing(
      sim_id = 1:n_sim,
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
      sim_id = 1:n_sim,
      all_data_grouped %>% select(dataset_label, x_cat),
      priors = factor(do_these_priors) %>% fct_inorder(), 
      value = NA_real_) %>%
    pivot_wider(id_cols = c(sim_id, dataset_label, priors),
                names_from = x_cat, 
                names_prefix = "alpha", 
                values_from = value,
                values_fill = NA_real_);
  
  xi_nans = 
    crossing(
      sim_id = 1:n_sim,
      all_data_grouped %>% select(dataset_label, x_cat),
      priors = factor(do_these_priors) %>% fct_inorder(), 
      value = NA_real_) %>%
    pivot_wider(id_cols = c(sim_id, dataset_label, priors),
                names_from = x_cat, 
                names_prefix = "xi", 
                values_from = value,
                values_fill = NA_real_);
  
  
  summarized_models = 
    crossing(
      all_data_grouped,
      sim_id = 1:n_sim,
      priors = factor(do_these_priors) %>% fct_inorder(), 
      #
      alphaj = NA_real_,
      #
      xij = NA_real_) %>%
    arrange(sim_id, dataset_label, priors, x)
  
  for(i in 1:n_sim) {
    
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
      curr_dir_shape = dir_shapes[[j]]
      
      for(curr_prior in do_these_priors) {
        
        cat("\n######################################\n");
        cat("#", curr_prior, "prior :: sim_id =", i, ":: dataset =", j,  "\n");
        cat("######################################\n\n");
        
        
        if(curr_prior %in% names(hs_stan_filenames)) {
          
          # horseshoe
          stan_args = 
            list(local_dof_stan = 1, 
                 global_dof_stan = 1,
                 alpha_scale_stan = curr_hs_scale[[curr_prior]],
                 slab_precision_stan = 1);
          
          stan_template = hs_stan_templates[[curr_prior]];
        } else if(curr_prior %in% names(ga_stan_filenames)) {
          
          # gamma
          stan_args = 
            list(
              alpha_shape_stan = curr_ga_shape[[curr_prior]],
              tiny_positive_stan = curr_ga_lower_bound[[curr_prior]]);
          
          stan_template = ga_stan_templates[[curr_prior]];
          
        } else {
          
          # dirichlet
          stan_args = 
            list(
              alpha_shape_stan = curr_dir_shape[[curr_prior]]);
          
          stan_template = dir_stan_templates[[curr_prior]];
          
        } 
        
        
        curr_row_performance = 
          with(summarized_performance, 
               which(dataset_label == j &
                       priors == curr_prior & 
                       sim_id == i));
        
        curr_row_nans = 
          with(alpha_nans, 
               which(dataset_label == j &
                       priors == curr_prior & 
                       sim_id == i));
        
        stopifnot(length(curr_row_performance) == 1)
        
        curr_rows_models <- 
          with(summarized_models,
               which(dataset_label == j &
                       priors == curr_prior & 
                       sim_id == i));
        
        if(curr_prior %in% names(hs_stan_filenames)) {
          
          summarized_performance[curr_row_performance, "tuning_param_val"] = 
            curr_hs_scale[[curr_prior]];
          summarized_performance[curr_row_performance, "lower_truncation"] = 
            0;
        } else if(curr_prior %in% names(ga_stan_filenames)) {
          summarized_performance[curr_row_performance, "tuning_param_val"] = 
            curr_ga_shape[[curr_prior]];
          summarized_performance[curr_row_performance, "lower_truncation"] = 
            curr_ga_lower_bound[[curr_prior]];
        } else {
          summarized_performance[curr_row_performance, "tuning_param_val"] = 
            curr_dir_shape[[curr_prior]];
          summarized_performance[curr_row_performance, "lower_truncation"] = 
            0
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
                                     stan_seed = stan_seeds[i],
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
                               priors == curr_prior & 
                               sim_id == i, curr_fit$alpha_number_nan[1:(curr_n_cats)], .)) 
        
        xi_nans <- 
          xi_nans %>%
          mutate_at(vars(paste0("xi",1:(curr_n_cats))),
                    ~ ifelse(dataset_label == j & 
                               priors == curr_prior & 
                               sim_id == i, curr_fit$xi_number_nan[1:(curr_n_cats)], .)) 
        
        summarized_models[curr_rows_models, "xij"] = 
          apply(curr_fit$all_draws$xi, 2, median)
        
        summarized_models[curr_rows_models, "alphaj"] = 
          apply(curr_fit$all_draws$alpha[, 1:curr_n_cats, drop = F], 2, median)
        
        
      } 
      
      write_csv(filter(summarized_performance, 
                       dataset_label <= j, 
                       sim_id <= i),
                path = "out/exemplar2_fixed_data_performance.csv",
                append = FALSE);
      
      write_csv(filter(summarized_models,
                       dataset_label <= j, 
                       sim_id <= i),
                path = "out/exemplar2_fixed_data_models.csv",
                append = FALSE);
      
      write_csv(filter(xi_nans, 
                       dataset_label <= j, 
                       sim_id <= i),
                path = "out/exemplar2_xi_nans.csv",
                append = FALSE);
      
      write_csv(filter(alpha_nans, 
                       dataset_label <= j, 
                       sim_id <= i),
                path = "out/exemplar2_alpha_nans.csv",
                append = FALSE);
      
      rm(stan_template, stan_args, curr_prior);
      rm(curr_row_performance, curr_row_nans, curr_rows_models);
      
      
    }
    rm(curr_ga_shape, curr_hs_scale, curr_ga_lower_bound);
    rm(data_grouped, curr_n_cats, curr_n_per_cat)
    
  }
  
} else {
  
  summarized_performance = 
    read_csv("out/exemplar2_fixed_data_performance.csv")
  
  summarized_models = 
    read_csv("out/exemplar2_fixed_data_models.csv")
  
  xi_nans = 
    read_csv("out/exemplar2_xi_nans.csv")
  
  alpha_nans =
    read_csv("out/exemplar2_alpha_nans.csv")
  
  
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
    filter(priors != "dirichlet1") %>%
    mutate(dataset_pretty_label = 
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
             ))
  
  
  ggplot(summarized_models) + 
    geom_boxplot(aes(x = factor(x), 
                     y = xij),
                 lwd = 0.25,
                 coef = 1e6, 
                 fill = grey(0.2),
                 color = grey(0.2)) + 
    geom_point(aes(x = factor(x), 
                   y = I(y/n)),
               size = 2.5,
               color = "sienna2",
               shape = "diamond") + 
    facet_grid(priors_pretty_label ~ dataset_pretty_label, 
               scales = "free_x",
               labeller = label_parsed) + 
    scale_x_discrete(name = "j") +
    scale_y_continuous(name = expression(xi[j])) +
    theme(text = element_text(size = 18));
  ggsave(filename = "numerical_xi.pdf", height = 12, width = 9)
  
  ggplot(filter(summarized_models, dataset_label == 1)) + 
    geom_boxplot(aes(x = factor(x), 
                     y = xij),
                 lwd = 0.25,
                 coef = 1e6, 
                 fill = grey(0.2),
                 color = grey(0.2)) + 
    geom_point(aes(x = factor(x), 
                   y = I(y/n)),
               size = 2.5,
               color = "sienna2",
               shape = "diamond") + 
    facet_wrap(~ priors_pretty_label, 
               labeller = label_parsed, 
               scales = "fixed") + 
    scale_x_discrete(name = "j") +
    scale_y_continuous(name = expression(xi[j])) +
    theme(text = element_text(size = 18));
  ggsave(filename = "numerical_xi_reduced.pdf", height = 7, width = 9)
  
  
  ggplot(summarized_models) + 
    geom_boxplot(aes(x = factor(x), 
                     y = log10(alphaj)),
                 lwd = 0.25,
                 coef = 1e6,
                 fill = grey(0.2),
                 color = grey(0.2)) + 
    facet_grid(priors_pretty_label ~ dataset_pretty_label, 
               labeller = label_parsed, 
               scales = "free") + 
    scale_x_discrete(name = "j") +
    scale_y_continuous(name = expression(log(alpha[j]))) +
    theme(text = element_text(size = 18));
  ggsave(filename = "numerical_alpha.pdf", height = 12, width = 9)
  
  ggplot(filter(summarized_models, dataset_label == 1)) + 
    geom_boxplot(aes(x = factor(x), 
                     y = log10(alphaj)),
                 lwd = 0.25,
                 coef = 1e6,
                 fill = grey(0.2),
                 color = grey(0.2)) + 
    facet_wrap(~ priors_pretty_label, 
               labeller = label_parsed, 
               scales = "fixed") + 
    scale_x_discrete(name = "j") +
    scale_y_continuous(name = expression(log(alpha[j]))) +
    theme(text = element_text(size = 18));
  ggsave(filename = "numerical_alpha_reduced.pdf", height = 7, width = 9)
  
  
  
}

