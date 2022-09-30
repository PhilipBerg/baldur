data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  int C;            // Number of comparisons to perform
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
  int c[C, 2];      // Comparisons to perform
  real alpha;       // alpha prior for gamma
  real beta_gamma;  // beta prior for gamma
  vector[K] xbar;   // prior mu
  vector[N] u;      // uncertainty
}
parameters {
  vector[K] mu;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
  real y_diff[C];       // difference in coefficients
  vector[K] eta;        // Error in mean
}
transformed parameters{
  vector[C] mu_diff;
  mu_diff = mu[c[, 1]] - mu[c[, 2]];
}
model {
  sigma  ~ gamma(alpha, beta_gamma);
  eta    ~ normal(0, 1);
  y      ~ normal(x * mu, sigma*u);
  mu     ~ normal(xbar + sigma*eta, sigma);
  y_diff ~ normal(mu_diff, sigma);
}
