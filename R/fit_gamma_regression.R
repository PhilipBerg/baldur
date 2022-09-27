#' Function for Fitting the Mean-Variance Gamma Regression Models
#'
#' `fit_gamma_regression` returns a `glm` object containing the
#'    gamma regression for the mean-variance trend.
#'
#' @param data a `data.frame` to generate the mean-variance trends for. It
#'     should contain columns with conditions named as the column names in
#'     `design` (presumably with some suffix).
#' @param formula a formula describing the model
#'
#' @param ... only for compatibility with other functions
#'
#' @return `fit_gamma_regression` returns a glm object
#' @export
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important for calculate_mean_sd_trends
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log-transform the data
#' yeast_norm <- psrn(yeast, "identifier") %>%
#'   # Add row means and variances
#'   calculate_mean_sd_trends(design)
#'
#' # Fit gamma regression model for the mean-variance trends
#' gamma_model <- fit_gamma_regression(yeast_norm, sd ~ mean)
fit_gamma_regression <- function(data, formula = sd ~ mean + c, ...) {
    stats::glm(formula, stats::Gamma(log), data)
}
