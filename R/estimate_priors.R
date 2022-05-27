#' Estimate Gamma priors for sigma
#'
#' @description Estimates the priors for the Bayesian model.
#' @param data A `tibble` or `data.frame` to add gamma priors to
#' @param design_matrix A design matrix for the data (see example)
#' @param formula Formula for the Gamma regression
#'
#' @return A `tibble` or `data.frame` with the alpha,beta priors
#' @export
#'
#' @examples 'lorem'
estimate_gamma_priors <- function(data, design_matrix, formula = sd ~ mean + c){
  gamma_reg <- glm(formula, Gamma(log), data)
  data %>%
    dplyr::mutate(
      alpha = (1/summary(gamma_reg)$dispersion),
      beta = estimate_beta(gamma_reg, mean, c, alpha)
    )
}
