data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of conditions
  int C;            // number of comparisons to perform
  matrix[N, K] x;   // design matrix
  vector[N] y;      // data
  int c[C, 2];      // Comparisons to perform
  real alpha;       // alpha prior for gamma
  real beta_gamma;  // beta prior for gamma
  vector[N] u;      // uncertainty
}
parameters {
  vector[K] mu;         // coefficients for predictors
  real<lower=0> sigma;  // error scale
  real y_diff[C];       // difference in coefficients
  vector[K] eta;        // Error in mean
  vector[K] mu_not;     // prior mu
}
transformed parameters{
  vector[C] mu_diff;    // differences in means
  mu_diff = mu[c[, 1]] - mu[c[, 2]];
}
model {
  sigma  ~ gamma(alpha, beta_gamma);        // variance
  eta    ~ normal(0, 1);                    // NCP auxilary variable
  mu     ~ normal(mu_not + sigma*eta, sigma); // mean
  y      ~ normal(x * mu, sigma*u);         // data model
  y_diff ~ normal(mu_diff, sigma);          // difference statistic
}
