# Save this file as `R/lm_stan.R`

#' Bayesian linear regression with Stan
#'
#' @export
#' @param x Numeric vector of input values.
#' @param y Numeric vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
lm_stan <- function(...) {
  standata <- list(N = 6, K = 2, C = 1, x = structure(c(1, 1, 1, 0, 0, 0, 0,
                                                        0, 0, 1, 1, 1), .Dim = c(6L, 2L), .Dimnames = list(c("1", "2",
                                                                                                             "3", "4", "5", "6"), c("x1", "x2")), assign = c(1L, 1L), contrasts = list(
                                                                                                               `factor(rep(c("x1", "x2"), each = 3))` = "contr.treatment")),
                   y = c(1.71731160937843, 2.30339453276428, 2.03280175657505,
                         1.53793954739669, 1.7368787487597, 1.9988433599649), c = structure(1:2, .Dim = 1:2),
                   alpha = 2.58401410965892, beta_gamma = 14.6326506797514,
                   xbar = c(x1 = 2.01783596623925, x2 = 1.7578872187071), u = c(x1_1 = 0.179056360610137,
                                                                                x1_2 = 0.170729962345089, x1_3 = 0.174524998982475, x2_1 = 0.181684935250435,
                                                                                x2_2 = 0.178771928084969, x2_3 = 0.175007181439148))
  out <- rstan::sampling(stanmodels$weighted_decision, data = standata, ...)
  return(out)
}
