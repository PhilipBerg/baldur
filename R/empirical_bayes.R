#' Baldur's empirical Bayes Prior For The Mean In Conditions
#'
#' @description Here we assume that the sample mean of each condition is an estimator for the center of the mean prior.
#' In addition, it assumes that the confidence in the prior is proportional to the variance of the peptide.
#'
#'
#' @section Code:
#' The `Stan` code for the model is given by:
#' ```stan
#' data {
#' int<lower=0> N;     // number of data items
#' int<lower=0> K;     // number of conditions
#' int C;              // number of comparisons to perform
#' matrix[N, K] x;     // design matrix
#' vector[N] y;        // data
#' matrix[K, C] c;     // contrast matrix
#' real alpha;         // alpha prior for gamma
#' real beta;          // beta prior for gamma
#' vector[N] u;        // uncertainty
#' vector[K] mu_not;   // prior mu
#' }
#' transformed data{
#'   vector[K] n_k;      // per condition measurements
#'   row_vector[C] n_c;  // per comparison measurements
#'   matrix[K, C] abs_c; // abs of C for n_c calculation
#'   for (i in 1:K) {
#'     n_k[i] = 1/sum(x[,i]);
#'     for (j in 1:C) {
#'       abs_c[i, j] = abs(c[i, j]);
#'     }
#'   }
#'   n_c = n_k' * abs_c;
#'   n_c = sqrt(n_c);
#' }
#' parameters {
#'   vector[K] mu;           // mean vector
#'   real<lower=0> sigma;    // error scale
#'   real y_diff[C];         // difference in coefficients
#'   vector[K] eta;          // Error in mean
#'   vector[K] prior_mu_not; // Estimation error
#' }
#' transformed parameters{
#'   row_vector[C] mu_diff = mu' * c;        // differences in means
#' vector[K] sigma_mu_not = 2*sigma * n_k; // variance of ybars
#' vector[C] sigma_lfc = sigma * n_c';     // variance of y_diff
#' }
#' model {
#'   sigma        ~ gamma(alpha, beta);                        // variance
#'   eta          ~ normal(0, 1);                              // NCP auxilary variable
#'   prior_mu_not ~ normal(mu_not, sigma_mu_not);              // Prior mean
#'   mu           ~ normal(prior_mu_not + sigma * eta, sigma); // mean
#'   y            ~ normal(x * mu, sigma * u);                 // data model
#'   y_diff       ~ normal(mu_diff, sigma_lfc);                // difference statistic
#' }
#' ```
#'
#' @export
#' @name empirical_bayes
NULL
