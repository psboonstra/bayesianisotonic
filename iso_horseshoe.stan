// Regularized Horseshoe Prior (stable version)
data {
  int<lower = 1> n_groups_stan; // number of unique covariate patterns
  int<lower = 0> n_per_group_stan[n_groups_stan]; // number observations per pattern
  int<lower = 0, upper = max(n_per_group_stan)> y_stan[n_groups_stan]; //
  real<lower = 0> local_dof_stan; // dof of pi(lambda), = 1
  real<lower = 0> global_dof_stan; // dof of pi(tau), = 1
  real<lower = 0> alpha_scale_stan; // c, Equation (8)2 
  int<lower = 0,upper = 1> only_prior_stan;//if 1, ignore the model and data and generate from the prior only
}
transformed data {
  real small_number = machine_precision()^0.5;
  real alpha_scale_stan_sq = square(alpha_scale_stan);
}
parameters {
  vector<lower = 0.0>[n_groups_stan+1] alpha_raw;
  // tau is decomposed into chi-square and inverse-gamma portions
  real<lower = 0.0> tau_glob_base_sq;
  real<lower = 0.0> tau_glob_scale_sq;
  // lambdas are decomposed into chi-square and inverse-gamma portions
  vector<lower = 0.0>[n_groups_stan+1] lambda_base_sq;//
  vector<lower = 0.0>[n_groups_stan+1] lambda_scale_sq;
}
transformed parameters {
  vector<lower = 0.0,upper = 1.0>[n_groups_stan] xi; // 
  vector<lower = 0.0,upper = 1.0>[n_groups_stan+1] theta;
  vector<lower = 0.0>[n_groups_stan+1] alpha;
  real<lower = 0.0> tau_glob_sq;//tau^2
  vector<lower = 0.0>[n_groups_stan+1] lambda_sq;//lambda^2
  vector<lower = 0.0, upper = 1.0>[n_groups_stan+1] normalized_alpha;
  tau_glob_sq = tau_glob_base_sq * tau_glob_scale_sq;
  lambda_sq = lambda_base_sq .* lambda_scale_sq;
  for(i in 1:(n_groups_stan+1)) {
    theta[i] = 1.0 / sqrt(1.0 + (1.0 / (alpha_scale_stan_sq * tau_glob_sq * lambda_sq[i])));
  }
  alpha = (theta .* alpha_raw);
  normalized_alpha = alpha / sum(alpha);
  xi[1] = normalized_alpha[1];
  if(n_groups_stan > 1) {
    for(i in 2:n_groups_stan) {
      xi[i] = xi[i-1] + normalized_alpha[i];
    }
  }
}
model {
  alpha_raw ~ normal(0.0, 1.0);
  tau_glob_base_sq ~ chi_square(1.0);
  tau_glob_scale_sq ~ inv_gamma(global_dof_stan/2.0, global_dof_stan/2.0);
  lambda_base_sq ~ chi_square(1.0);
  lambda_scale_sq ~ inv_gamma(local_dof_stan/2.0, local_dof_stan/2.0);
  if(only_prior_stan == 0) {
    y_stan ~ binomial(n_per_group_stan, xi);
  }
}
