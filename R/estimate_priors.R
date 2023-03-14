#' Estimate Gamma priors for sigma
#'
#' @description Estimates the priors for the Bayesian data and decision model.
#'   `estimate_gamma_priors` is a wrapper that adds new columns to the data (one
#'   for alpha and one for betas).
#' @param data A `tibble` or `data.frame` to add gamma priors to
#' @param reg A `glm` Gamma regression or `lgmr` object
#'
#' @return `estimate_gamma_priors` returns a `tibble` or `data.frame` with the alpha,beta priors estimate.
#' @export
#'
#' @name estimate_gamma_priors
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
#'     estimate_gamma_priors(gam_reg)
#' # Can also explicitly estimate the beta parameters
#' estimate_beta(gam_reg, yeast_norm$mean, 1/summary(gamma_reg)$dispersion)
estimate_gamma_priors <- function(data, reg){
  if (!"mean" %in% names(data)) {
    stop('Mean column missnig\nDid you forget to calculate the M-V trend?\n Try running calculate_mean_sd_trends')
  }
  data %>%
    dplyr::mutate(
      alpha = 1/summary(reg)$dispersion,
      beta = estimate_beta(reg, mean, alpha)
    )
}

#' Estimate the beta parameter of the gamma distribution of sigma
#'
#' @param mean  The mean value of the peptide
#' @param alpha The alpha parameter of the peptide
#'
#' @rdname estimate_gamma_priors
#'
#' @return `estimate_gamma_priors` returns estimates of the beta parameter(s)
#' @export
estimate_beta <- function(mean, reg, alpha){
  alpha / stats::predict.glm(reg, newdata = data.frame(mean = mean), type = 'response')
}
