data {
  int<lower=0> N;     // number of data items
  int<lower=0> K;     // number of conditions
  int C;              // number of comparisons to perform
  matrix[N, K] x;     // design matrix
  vector[N] y;        // data
  matrix[C, K] c;     // contrast matrix
  real alpha;         // alpha prior for gamma
  real beta_gamma;    // beta prior for gamma
  vector[N] u;        // uncertainty
  vector[K] mu_not;   // prior mu
}

transformed data{
  vector[K] n_k;      // per condition measurements
  vector[C] n_c;      // per comparison measurements
  matrix[C, K] abs_c; // abs of C for n_c calculation
  for (i in 1:K) {
    n_k[i] = 1/sum(x[,i]);
    for (j in 1:C) {
      abs_c[j, i] = abs(c[j, i]);
    }
  }
  n_c = abs_c * n_k;
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
  vector[C] mu_diff;      // differences in means
  vector[K] sigma_mu_not; // variance of ybars
  vector[C] sigma_lfc;    // variance of y_diff
  mu_diff = c*mu;
  sigma_mu_not = 2*sigma * n_k;
  sigma_lfc = sigma * n_c;
}

model {
  sigma        ~ gamma(alpha, beta_gamma);                    // variance
  eta          ~ normal(0, 1);                                // NCP auxilary variable
  prior_mu_not ~ normal(mu_not, sigma_mu_not);                // Prior mean
  mu           ~ normal(prior_mu_not + sigma*eta, sigma);     // mean
  y            ~ normal(x * mu, sigma*u);                     // data model
  y_diff       ~ normal(mu_diff, sigma_lfc);                  // difference statistic
}
