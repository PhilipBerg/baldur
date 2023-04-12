data {
  int<lower=0> N;     // number of data items
  int<lower=0> K;     // number of conditions
  int C;              // number of comparisons to perform
  matrix[N, K] x;     // design matrix
  vector[N] y;        // data
  matrix[K, C] c;     // contrast matrix
  real alpha;         // alpha prior for gamma
  real beta;          // beta prior for gamma
  vector[N] u;        // uncertainty
}

transformed data{
  vector[K] n_k;                // per condition measurements
  row_vector[C] n_c;            // per comparison measurements
  matrix[K, C] abs_c = fabs(c); // abs of C for n_c calculation
  for (i in 1:K) {
    n_k[i] = 1/sum(x[,i]);
  }
  n_c = n_k' * abs_c;
  n_c = sqrt(n_c);
}

parameters {
  vector[K] mu;           // coefficients for predictors
  real<lower=0> sigma;    // error scale
  real y_diff[C];         // difference in coefficients
  vector[K] eta;          // Error in mean
  vector[K] prior_mu_not; // Estimation error
}

transformed parameters{
  row_vector[C] mu_diff = mu' * c;        // differences in means
  vector[C] sigma_lfc = sigma * n_c';     // variance of y_diff
}

model {
  sigma        ~ gamma(alpha, beta);                      // variance
  eta          ~ normal(0, 1);                            // NCP auxilary variable
  prior_mu_not ~ normal(0, 10);                          // prior mean
  mu           ~ normal(prior_mu_not + sigma*eta, sigma); // mean
  y            ~ normal(x * mu, sigma*u);                 // data model
  y_diff       ~ normal(mu_diff, sigma_lfc);              // difference statistic
}
