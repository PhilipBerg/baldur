#' Estimate measurement uncertainty
#' @description Estimates the measurement uncertainty for each data point using
#'   a Gamma regression.
#'
#' @param reg A `glm` gamma regression or `lgmr` object
#' @param data A `tibble` or `data.frame`
#' @param id_col A character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#' @param design_matrix Cell mean design matrix for the data
#'
#' @return A matrix with the uncertainty
#' @export
#'
#' @importFrom stats predict.glm
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr select
#' @importFrom dplyr matches
#' @importFrom magrittr set_rownames
#'
#' @name estimate_uncertainty
#'
#' @examples
#' # Setup model matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast_norm <- yeast %>%
#'   # Remove missing data
#'   tidyr::drop_na() %>%
#'   # Normalize data
#'   psrn("identifier") %>%
#'   # Add mean-variance trends
#'   calculate_mean_sd_trends(design)
#' # Fit the gamma regression
#' gam <- fit_gamma_regression(yeast_norm, sd ~ mean)
#' # Estimate each data point's uncertainty
#' estimate_uncertainty(gam, yeast_norm, 'identifier', design)

estimate_uncertainty <- function(reg, data, id_col, design_matrix){
  UseMethod("estimate_uncertainty")
}

#' @rdname estimate_uncertainty
#'
#' @export
estimate_uncertainty.glm <- function(reg, data, id_col, design_matrix){
  check_id_col(id_col, colnames(data))
  pred <- ~ stats::predict.glm(
    reg,
    newdata = data.frame(mean = .),
    type = "response"
  )
  condi_regex <- colnames(design_matrix) %>%
    paste0(collapse = '|')

  data %>%
    dplyr::mutate(
      dplyr::across(dplyr::matches(condi_regex),
                    pred
      )
    ) %>%
    dplyr::select(dplyr::matches(condi_regex)) %>%
    as.matrix() %>%
    magrittr::set_rownames(data[[id_col]])
}

#' @rdname estimate_uncertainty
#'
#' @export
estimate_uncertainty.lgmr <- function(reg, data, id_col, design_matrix){
  check_id_col(id_col, colnames(data))
  pars <- coef.lgmr(reg, simplify = TRUE, pars = c('coef', 'theta'))
  condi_regex <- colnames(design_matrix) %>%
    paste0(collapse = '|')

  theta <- sort_theta(pars$theta, data[[id_col]])

  mu_inputs <- mu_std_inputs(data)

  data %>%
    dplyr::mutate(
      dplyr::across(dplyr::matches(condi_regex),
                    ~ mu_fun(theta, pars$coef, ., mu_inputs[1], mu_inputs[2])
      )
    ) %>%
    dplyr::select(dplyr::matches(condi_regex)) %>%
    as.matrix() %>%
    magrittr::set_rownames(data[[id_col]])
}

sort_theta <- function(theta, reference) {
  names(theta) <- stringr::str_remove(names(theta), 'theta_')
  theta[match(reference, names(theta))]
}
