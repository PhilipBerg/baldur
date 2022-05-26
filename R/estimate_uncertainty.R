#' Estimate measurement uncertainty
#' @description Estimates the measurement uncertainty for each data point.
#' @param data A `tibble` or `data.frame`
#' @param formula Formula for the gamma regression
#' @param identifier Id column
#'
#' @return A matrix with the uncertainty
#' @export
#'
#' @examples 'lorem'
estimate_uncertainty <- function(data, identifier, formula = sd ~ mean + c){
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
  data %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric),
                    !!pred
      )
    ) %>%
    dplyr::select(-c(dplyr::all_of(identifier), -where(is.numeric), c, mean, sd)) %>%
    as.matrix() %>%
    magrittr::set_rownames(data[[identifier]])
}
