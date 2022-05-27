#' Estimate measurement uncertainty
#' @description Estimates the measurement uncertainty for each data point using a Gamma regression.
#' @param data A `tibble` or `data.frame`
#' @param identifier Id column
#' @param design_matrix Cell mean design matrix for the data
#' @param formula Formula for the gamma regression
#'
#' @return A matrix with the uncertainty
#' @export
#'
#' @examples 'lorem'
estimate_uncertainty <- function(data, identifier, design_matrix, formula = sd ~ mean + c){
  gamma_reg <- glm(formula, Gamma(log), data)
  if(!is.null(data$c)){
    pred <- ~predict.glm(
      gamma_reg,
      newdata = data.frame(mean = ., c = c),
      type = "response"
    )
  }else{
    pred <- ~predict.glm(
      gamma_reg,
      newdata = data.frame(mean = .),
      type = "response"
    )
  }
  condi_regex <- colnames(design_matrix) %>%
    paste0(collapse = '|')
  data %>%
    dplyr::mutate(
      dplyr::across(dplyr::where(is.numeric),
                    !!pred
      )
    ) %>%
    dplyr::select(dplyr::matches(condi_regex)) %>%
    as.matrix() %>%
    magrittr::set_rownames(data[[identifier]])
}
