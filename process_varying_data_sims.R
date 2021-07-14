
# Read in, process ----

library(tidyverse);library(knitr);library(kableExtra);
source("varying_data_functions.R");
source("functions.R");

col_types <- cols(
  priors = col_character(),
  sim_id = col_integer(),
  array_id = col_integer(), 
  scenario_id = col_integer(),
  true_prob_curve_id = col_integer(), 
  predictor_dist_id = col_integer(), 
  n_training = col_integer(), 
  rmse = col_double(),
  bias = col_double(),    
  loglik = col_double(), 
  kl_div = col_double(),
  loglik_intercept = col_double(),
  coverage_50 = col_double()
);

col_types_bayesian <- cols(
  priors = col_character(),
  sim_id = col_integer(),
  array_id = col_integer(),
  scenario_id = col_integer(),
  true_prob_curve_id = col_integer(),
  predictor_dist_id = col_integer(),
  n_training = col_integer(),
  number_divergences = col_integer(),
  any_nan = col_integer()
);

raw_all_performance =
  raw_all_bayesian_performance = tibble();
for(i in 1:(800)) {
  foo <- try(read_csv(paste0("out/job",i,"_performance.csv"), col_types = col_types));
  if(!"try-error" %in% class(foo)) {
    raw_all_performance = 
      bind_rows(raw_all_performance, foo);
  } else {
    cat("sim ", i, ", not found\n");
  }
  foo <- try(read_csv(paste0("out/job",i,"_bayesian_performance.csv"), col_types = col_types_bayesian));
  if(!"try-error" %in% class(foo)) {
    raw_all_bayesian_performance = 
      bind_rows(raw_all_bayesian_performance, foo);
  } 
}



all_performance <- 
  full_join(raw_all_performance,
            select(raw_all_bayesian_performance,
                   priors:n_training,tuning_param_val,
                   number_divergences:chain_relative_diff, 
                   kl_div) %>%
              rename(bayes_kl_div = kl_div),
            by = c("priors", "sim_id", "array_id",
                   "scenario_id", "true_prob_curve_id", "predictor_dist_id",
                   "n_training")) %>%
  mutate(scenario_id = factor(scenario_id),
         true_prob_curve_id = factor(true_prob_curve_id),
         predictor_dist_id = factor(predictor_dist_id))


priors_pretty = 
  c(nonbayes1 = "Isoreg",
    nonbayes2 = "IsoregMod",
    horseshoe = "HS", 
    gamma1 = "GA$_1$", 
    gamma2 = "GA$_2$", 
    gamma3 = "GA$_3$",
    gamma4 = "GA$_4$") %>% 
  factor() %>%
  fct_inorder();


#number sims
all_performance %>%
  group_by(true_prob_curve_id, predictor_dist_id, n_training,priors) %>%
  count() %>%
  group_by(true_prob_curve_id, predictor_dist_id, n_training) %>%
  summarize(number_max_sims = max(n),
            number_min_sims = min(n)) %>%
  arrange(true_prob_curve_id, predictor_dist_id, n_training) %>%
  as.data.frame();


# Table 2 Pointwise KL divergence ----

## + setup ----
priors_to_include = 
  c("nonbayes1",
    "horseshoe",
    "gamma1",
    "gamma2",
    "gamma3",
    "gamma4");

curr_data <- 
  all_performance %>%
  filter(priors %in% priors_to_include) %>%
  mutate(priors = 
           factor(priors, 
                  levels = priors_to_include,
                  labels = priors_pretty[priors_to_include], 
                  ordered = T));


dim(curr_data);

table2_colnames = c("Curve", 
                    "$n$",
                    "Metric",
                    levels(curr_data$priors));

table2_header = c("Design factors" = 2, 
                  " " = 1,
                  "Methods" = length(levels(curr_data$priors)));


table2_tall <- 
  curr_data %>%
  group_by(true_prob_curve_id, 
           n_training, 
           priors, 
  ) %>%
  summarize(mean_kl_div = mean(1e3 * kl_div),
            mean_rmse = mean(1e3 * rmse)) %>%
  group_by(true_prob_curve_id,  n_training) %>%
  mutate(best_kl_div2 = ifelse(mean_kl_div <= 1.10 * min(mean_kl_div),"\\textbf{", ""),
         best_kl_div3 = ifelse(mean_kl_div <= 1.10 * min(mean_kl_div),"}", ""),
         best_rmse2 = ifelse(mean_rmse <= 1.10 * min(mean_rmse),"\\textbf{", ""),
         best_rmse3 = ifelse(mean_rmse <= 1.10 * min(mean_rmse),"}", "")) %>%
  ungroup() %>%
  arrange(true_prob_curve_id, n_training, priors) %>%
  mutate(
    mean_kl_div = formatC(mean_kl_div, format = "f", digits = 0),
    pretty_kl_div = paste0(best_kl_div2,mean_kl_div, best_kl_div3),
    mean_rmse = formatC(mean_rmse, format = "f", digits = 0),
    pretty_rmse = paste0(best_rmse2,mean_rmse, best_rmse3)) %>%
  select(true_prob_curve_id, 
         n_training, 
         priors, 
         pretty_kl_div,
         pretty_rmse);

table2 <- 
  bind_rows(
    bind_cols(
      metric = "KL Divergence",
      table2_tall %>% 
        pivot_wider(
          id_cols = true_prob_curve_id:n_training, 
          names_from = priors,          
          values_from = pretty_kl_div)),
    bind_cols(
      metric = "RMSE",
      table2_tall %>% 
        pivot_wider(
          id_cols = true_prob_curve_id:n_training, 
          names_from = priors,          
          values_from = pretty_rmse))) %>%
  select(true_prob_curve_id,
         n_training, 
         everything())

linesep_index <- rep("", nrow(table2));
linesep_index[4] = "\\addlinespace";
#linesep_index[!duplicated(table2[,c("true_prob_curve_id")], fromLast = T)] = "\\addlinespace";


table2 %>%
  knitr::kable(format = "latex",
               col.names = table2_colnames,
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = linesep_index) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(table2_header);


# Figs 1, 2: Sample curves ----
# setup ----

col_types <- cols(
  array_id = col_integer(), 
  scenario_id = col_integer(),
  true_prob_curve_id = col_integer(), 
  predictor_dist_id = col_integer(), 
  n_training = col_integer(), 
  sim_id = col_integer(),
  priors = col_character(),
  x = col_double(),
  true_prob= col_double(),     
  fitted_prob= col_double()
);

raw_all_models = tibble();
for(i in 1:(800)) {
  foo <- try(read_csv(paste0("out/job",i,"_models.csv"), col_types = col_types));
  if(!"try-error" %in% class(foo)) {
    raw_all_models = 
      bind_rows(raw_all_models, 
                bind_cols( job_id = rep(i, nrow(foo)), foo));
  } else {
    cat("sim ", i, ", not found\n");
  }
}
rm(foo);

priors_to_include = c(
  "nonbayes1",
  "horseshoe",
  "gamma1",
  #"gamma2",
  #"gamma3",
  "gamma4")

all_models <- 
  raw_all_models %>%
  filter(priors %in% priors_to_include) %>%
  arrange(n_training, predictor_dist_id, true_prob_curve_id, array_id, sim_id, x) %>%
  mutate(n_training = paste0("n = ", n_training) %>% factor() %>% fct_inorder(), 
         bias = fitted_prob - true_prob, 
         scenario_id = factor(scenario_id),
         true_prob_curve_id = paste0("Curve ", true_prob_curve_id) %>% factor() %>% fct_inorder(),
         predictor_dist_id = factor(predictor_dist_id));

subset_array_ids <-
  all_models %>%
  filter(predictor_dist_id == 1) %>%
  group_by(scenario_id) %>%
  sample_n(4) %>%
  ungroup() %>%
  select(array_id, sim_id) %>%
  arrange(array_id, sim_id)


all_models_random_subset <- 
  all_models %>%
  #right_join(subset_array_ids) %>%
  arrange(scenario_id, array_id, sim_id) %>%
  mutate(n_training = fct_drop(n_training)) %>%
  left_join(all_performance %>%
              select(array_id, sim_id, 
                     true_prob_curve_id, predictor_dist_id, n_training, priors, 
                     kl_div, number_divergences) %>%
              mutate(n_training = paste0("n = ", n_training) %>% factor() %>% fct_inorder())) %>%
  mutate(priors = 
           factor(priors, 
                  levels = priors_to_include,
                  labels = priors_pretty[priors_to_include], 
                  ordered = T),
         n_training = fct_drop(n_training));

ggplot() +
  geom_boxplot(data = all_models_random_subset,
               aes(x = factor(x), 
                   y = fitted_prob, 
                   color = priors,
                   fill = priors),
               position = "dodge",
               coef = 1e6) +
  #geom_line(data = all_models_randdom_subset,
  #          aes(x = factor(x), 
  #              y = fitted_prob, 
  #              color = priors,
  #              group = interaction(priors, array_id, sim_id)),
  #          alpha = 0.4, 
  #          size = 1.1) + 
  geom_point(data =
               filter(all_models,
                      predictor_dist_id == 1,
                      sim_id == 1) %>%
               arrange(true_prob_curve_id, array_id) %>%
               group_by(true_prob_curve_id, n_training) %>%
               filter(array_id == dplyr::first(array_id)),
             aes(x = factor(x), 
                 y = true_prob), 
             color = "black",
             alpha = 0.9,
             size = 1.5) +
  facet_grid(n_training ~ true_prob_curve_id, scales = "free_y") + 
  coord_cartesian(ylim = c(0,1)) + 
  scale_y_continuous(name = "Fitted / true probability") + 
  scale_x_discrete(name = "X", labels = NULL) + 
  scale_color_brewer(name = "Prior",
                     palette = "Dark2",
                     labels = c(
                       "Isoreg" = "Isoreg",
                       "HS" = "HS",
                       "GA$_1$" = expression(GA[1]),
                       "GA$_2$" = expression(GA[2]),
                       "GA$_3$" = expression(GA[3]),
                       "GA$_4$" = expression(GA[4]))) + 
  scale_fill_brewer(name = "Prior",
                    palette = "Dark2",
                    labels = c(
                      "Isoreg" = "Isoreg",
                      "HS" = "HS",
                      "GA$_1$" = expression(GA[1]),
                      "GA$_2$" = expression(GA[2]),
                      "GA$_3$" = expression(GA[3]),
                      "GA$_4$" = expression(GA[4]))) + 
  #guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.7)),
  #       fill = guide_legend(override.aes = list(alpha = 1, size = 1.7))) + 
  theme(legend.position = "top", 
        legend.direction = "horizontal", 
        legend.key.width = unit(1.5,"cm"),
        text = element_text(size = 20));
ggsave(filename = "../fig1.pdf", height = 8, width = 11);

rm(priors_to_include);

