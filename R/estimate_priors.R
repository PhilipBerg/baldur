#' Estimate Gamma hyperparameters for sigma
#'
#' @description Estimates the hyperparameters for the Bayesian data and decision
#'   model. `estimate_gamma_hyperparameters` is a wrapper that adds new columns
#'   to the data (one for alpha and one for betas). Note that for `lgmr`
#'   objects, the `estimate_beta` function assumes that the data is ordered as
#'   when the model was fitted. If this is not the case, theta's will be
#'   incorrectly matched with peptides---resulting in wrong estimates of beta
#'   parameters. On the other hand, `estimate_gamma_hyperparameters` will (if
#'   needed) temporarily sort the data as when fitted and the sort it back as it
#'   was input to the function.
#' @param reg A `glm` Gamma regression or a `lgmr` object
#' @param data A `tibble` or `data.frame` to add gamma priors to
#' @param mean  The mean value of the peptide
#' @param m The mean of the means
#' @param s The sd of the means
#' @param alpha The alpha parameter of the peptide
#' @param ... Currently not in use
#'
#' @return `estimate_gamma_hyperparameters` returns a `tibble` or `data.frame`
#'   with the alpha,beta hyperparameters estimates as new columns.
#'
#' @export
#'
#'
#' @name estimate_gamma_hyperparameters
#'
#' @importFrom dplyr mutate
#' @importFrom stats predict.glm
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
#'
#' # Fit gamma regression (could also have been a lgmr model)
#' gam_reg <- fit_gamma_regression(yeast_norm, sd ~ mean)
#'
#' # Estimate priors
#' gam_reg %>%
#'     estimate_gamma_hyperparameters(yeast_norm)
#'
#' # Can also explicitly estimate the beta parameters
#' # Note this is order sensitive.
#' estimate_beta(gam_reg, yeast_norm$mean, 1/summary(gam_reg)$dispersion)
estimate_gamma_hyperparameters <- function(reg, data){
  UseMethod("estimate_gamma_hyperparameters")
}

#' @rdname estimate_gamma_hyperparameters
#' @export
estimate_gamma_hyperparameters.glm <- function(reg, data){
  if (!"mean" %in% names(data)) {
    stop("Mean column missnig\nDid you forget to calculate the M-V trend?\n Try running calculate_mean_sd_trends")
  }
  data %>%
    dplyr::mutate(
      alpha = 1/summary(reg)$dispersion,
      beta = estimate_beta(reg, mean, alpha)
    )
}

#' @rdname estimate_gamma_hyperparameters
#' @export
estimate_gamma_hyperparameters.lgmr <- function(reg, data){

  mu_inputs <- mu_std_inputs(data)
  pars <- coef(reg, simplify = TRUE, pars = c("all"))

  ori_order <- data$mean
  data <- orderer(data, reg$data$mean)

  data %>%
    dplyr::mutate(
      alpha = pars$aux["alpha"],
      beta  = alpha / mu_fun(pars$theta, pars$coef, mean, mu_inputs[1], mu_inputs[2])
    ) %>%
    orderer(ori_order)
}

#' @rdname estimate_gamma_hyperparameters
#'
#' @return `estimate_beta` returns estimates of the beta parameter(s)
#' @export
estimate_beta <- function(reg, mean, alpha, m, s, ...){
  UseMethod("estimate_beta")
}

#' @rdname estimate_gamma_hyperparameters
#'
#' @export
estimate_beta.glm <- function(reg, mean, alpha, m, s, ...){
  alpha / stats::predict.glm(reg, newdata = data.frame(mean = mean), type = 'response')
}

#' @rdname estimate_gamma_hyperparameters
#'
#' @export
estimate_beta.lgmr <- function(reg, mean, m, s, ...){
  pars <- coef(reg, simplify = TRUE, pars = c("all"))
  pars$aux["alpha"] / mu_fun(pars$theta, pars$coef, mean, m, s)
}

orderer <- function(data, order) {
  check_order <- any(data$mean != order)
  if (check_order) {
    data <- data[match(order, data$mean),]
  }
  return(data)
}

mu_std_inputs <- function(data) {

  if (!"mean" %in% names(data)) {
    stop("Mean column missnig\nDid you forget to calculate the M-V trend?\n Try running calculate_mean_sd_trends")
  }

  m <- mean(data$mean)
  s <- sd(data$s)

  c(m, s)
}
