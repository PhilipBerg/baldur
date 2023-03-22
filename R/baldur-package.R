#' The 'baldur' package.
#'
#' @description baldur is a Bayesian hierarchical model for statistical decision
#'   in proteomics data. It partitions the mean-variance trend into two parts
#'   using a probabilistic distance metric based on a gamma regression. It then
#'   uses a gamma regression (with or without partitioning) as en Empirical
#'   Bayes estimator for the prior on the variance. Further, it assumes that
#'   each measurement has an uncertainty (increased variance) associated with it
#'   that it also infers. Finally, it tries to estimate the posterior
#'   distribution (by Hamiltonian Monte Carlo) for the differences in means for
#'   each peptide in the data. Once the posterior is inferred, it integrates the
#'   tails to estimate the probability of error from which a statistical
#'   decision can be made.
#'
#' @docType package
#' @name baldur-package
#' @aliases baldur
#' @useDynLib baldur, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom stats Gamma
#' @importFrom stats glm
#' @importFrom stats predict.glm
#' @importFrom stats residuals
#' @importFrom stats sd
#' @importFrom stats dgamma
#' @importFrom magrittr %>%
#' @importFrom magrittr set_rownames
#' @importFrom dplyr across
#' @importFrom dplyr mutate
#' @importFrom dplyr matches
#' @importFrom rlang dots_list
#' @references Berg and Popescu (2023). Baldur: a Bayesian hierarchical model for
#' label-free proteomics using mean-variance trend. Under review. Stan
#' Development Team (2022). RStan: the R interface to Stan. R package version
#' 2.21.5. https://mc-stan.org
#'
NULL
