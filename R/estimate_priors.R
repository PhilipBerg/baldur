#' Estimate Gamma priors for sigma
#'
#' @description Estimates the priors for the Bayesian model.
#' @param data A `tibble` or `data.frame` to add gamma priors to
#' @param design_matrix A design matrix for the data (see example)
#' @param gamma_reg A Gamma regression object
#'
#' @return A `tibble` or `data.frame` with the alpha,beta priors
#' @export
#'
#' @importFrom dplyr mutate
#'
#' @examples
#' # Setup model matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize data
#' yeast_norm <- yeast %>%
#'     psrn("identifier") %>%
#'     # Get mean-variance trends
#'     calculate_mean_sd_trends(design)
#' # Fit gamma regression
#' gam_reg <- fit_gamma_regression(yeast_norm, sd ~ mean)
#'
#' # Estimate priors
#' yeast_norm %>%
#'     estimate_gamma_priors(design, gam_reg)
estimate_gamma_priors <- function(data, design_matrix, gamma_reg){
  data %>%
    dplyr::mutate(
      alpha = 1/summary(gamma_reg)$dispersion,
      beta = estimate_beta(gamma_reg, mean, c, alpha)
    )
}
