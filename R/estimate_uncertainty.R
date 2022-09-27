#' Estimate measurement uncertainty
#' @description Estimates the measurement uncertainty for each data point using a Gamma regression.
#' @param data A `tibble` or `data.frame`
#' @param identifier Name of the peptide name/id column
#' @param design_matrix Cell mean design matrix for the data
#' @param gam_reg A gamma regression model
#'
#' @return A matrix with the uncertainty
#' @export
#'
#' @examples
#' #' # Setup model matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast_norm <- yeast %>%
#' # Remove missing data
#' tidyr::drop_na() %>%
#'   # Normalize data
#' psrn('identifier') %>%
#' # Add mean-variance trends
#' calculate_mean_sd_trends(design)
#' # Fit the gamma regression
#' gam <- fit_gamma_regression(yeast_norm, sd ~ mean)
#' # Estimate each data point's uncertainty
#' unc <- estimate_uncertainty(yeast_norm, 'identifier', design, gam)
estimate_uncertainty <- function(data, identifier, design_matrix, gam_reg){
  if('c' %in% names(data)){
    pred <- ~ stats::predict.glm(
      gam_reg,
      newdata = data.frame(mean = ., c = c),
      type = "response"
    )
  } else {
    pred <- ~ stats::predict.glm(
      gam_reg,
      newdata = data.frame(mean = .),
      type = "response"
    )
  }
  condi_regex <- colnames(design_matrix) %>%
    paste0(collapse = '|')
  data %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric),
                    !!pred
      )
    ) %>%
    dplyr::select(dplyr::matches(condi_regex)) %>%
    as.matrix() %>%
    magrittr::set_rownames(data[[identifier]])
}
