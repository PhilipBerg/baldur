data {
  int<lower=0> N;     // number of data items
  int<lower=0> K;     // number of conditions
  int C;              // number of comparisons to perform
  matrix[N, K] x;     // design matrix
  vector[N] y;        // data
  int c[C, 2];        // Comparisons to perform
  real alpha;         // alpha prior for gamma
  real beta_gamma;    // beta prior for gamma
  vector[N] u;        // uncertainty
}
transformed data{
  vector[K] n_k;     // per condition measurements
  vector[C] n_c;     // per comparison measurements
  for (i in 1:K) {
    n_k[i] = 1/sum(x[,i]);
  }
  n_c = sqrt((n_k[c[, 1]] + n_k[c[, 2]]));
}
parameters {
  vector[K] mu;           // coefficients for predictors
  real<lower=0> sigma;    // error scale
  real y_diff[C];         // difference in coefficients
  vector[K] eta;          // Error in mean
  vector[K] prior_mu_not; // Estimation error
}
transformed parameters{
  vector[C] mu_diff;      // differences in means
  vector[C] sigma_lfc;    // variance of ybars
  mu_diff = mu[c[, 1]] - mu[c[, 2]];
  sigma_lfc = sigma * n_c;
}
model {
  sigma        ~ gamma(alpha, beta_gamma);                    // variance
  eta          ~ normal(0, 1);                                // NCP auxilary variable
  prior_mu_not ~ normal(0, 100);                              // prior mean
  mu           ~ normal(prior_mu_not + sigma*eta, sigma);     // mean
  y            ~ normal(x * mu, sigma*u);                     // data model
  y_diff       ~ normal(mu_diff, sigma_lfc);                  // difference statistic
}
