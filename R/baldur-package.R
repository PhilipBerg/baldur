#' The 'baldur' package.
#'
#' @description `baldur` is a Bayesian hierarchical model for statistical decision
#'   in proteomics data. It models the mean-variance trend with the option of
#'   two different regression models, a gamma regression or a latent gamma
#'   mixture regression. It then the regression model as en Empirical Bayes
#'   estimator for the prior on the variance. Further, it assumes that
#'   each measurement has an uncertainty (increased variance) associated with it
#'   that it also infers. Finally, it tries to estimate the posterior
#'   distribution (by Hamiltonian Monte Carlo) for the differences in means for
#'   each peptide in the data. Once the posterior is inferred, it integrates the
#'   tails to estimate the probability of error from which a statistical
#'   decision can be made.
#'
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
#' @importFrom Rdpack reprompt
#' @references Berg George Popescu (2023) "Baldur: Bayesian Hierarchical Modeling for Label-Free Proteomics with Gamma Regressing Mean-Variance Trends"
#' Molecular & Cellular Proteomics: 2023-12. https://doi.org/10.1016/j.mcpro.2023.100658
#'
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version
#' 2.21.5. https://mc-stan.org
#'
## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
"_PACKAGE"
