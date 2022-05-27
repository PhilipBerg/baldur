#' Calculate the Mean-Variance trend
#'
#' @description Calculates the mean and sd of the rows.
#' @param data A `tibble` or `data.frame` to annotate with mean and sd
#' @param design_matrix A design matrix for the data (see example)
#'
#' @return A `tibble` or `data.frame` with the mean and sd vectors
#' @export
#'
#' @examples 'lorem'
calculate_mean_sd_trends <- function(data, design_matrix){
  conditions <- design_matrix %>%
    colnames() %>%
    paste0(collapse = '|')
  data %>%
    dplyr::mutate(
      mean = rowMeans(dplyr::across(dplyr::matches(conditions))),
      sd = apply(dplyr::across(dplyr::matches(conditions)), 1, sd)
    )
}
