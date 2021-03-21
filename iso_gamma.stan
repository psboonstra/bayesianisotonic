// Dirichlet distribution on each of the increments
data {
  int<lower = 1> n_groups_stan; // number of unique covariate patterns
  int<lower = 0> n_per_group_stan[n_groups_stan]; // number observations per pattern
  int<lower = 0, upper = max(n_per_group_stan)> y_stan[n_groups_stan]; // 
  real<lower = 0> alpha_shape_stan;//
  real<lower = 0> tiny_positive_stan;
  int<lower = 0,upper = 1> only_prior_stan;//if 1, ignore the model and data and generate from the prior only
}
transformed data {
  real scaled_tiny_positive_stan = ((n_groups_stan + 1) * tiny_positive_stan);
}
parameters {
  vector<lower = tiny_positive_stan>[n_groups_stan + 1] alpha; //
}
transformed parameters {
  vector<lower = 0.0,upper = 1.0>[n_groups_stan] xi; // 
  real<lower = scaled_tiny_positive_stan> sum_alpha;
  sum_alpha =  scaled_tiny_positive_stan + (1 + scaled_tiny_positive_stan) * sum(alpha);
  xi[1] = alpha[1] / sum_alpha;
  if(n_groups_stan > 1) {
    for(i in 2:n_groups_stan) {
      xi[i] = xi[i-1] + (alpha[i] / sum_alpha);
    }
  }
}
model {
  alpha ~ gamma(alpha_shape_stan, 1.0);
  if(only_prior_stan == 0) {
    y_stan ~ binomial(n_per_group_stan, xi);
  }
}
