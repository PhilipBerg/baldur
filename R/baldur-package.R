#' The 'baldur' package.
#'
#' @description baldur is a Bayesian heirachal model for statistical decision in proteomics data
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
#' @importFrom dplyr across
#' @references
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version 2.21.5. https://mc-stan.org
#'
NULL
