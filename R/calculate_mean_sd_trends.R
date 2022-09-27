#' Calculate the Mean-Variance trend
#'
#' @description Calculates the mean and sd of the rows.
#' @param data A `tibble` or `data.frame` to annotate with mean and sd
#' @param design_matrix A design matrix for the data (see example)
#'
#' @return A `tibble` or `data.frame` with the mean and sd vectors
#' @export
#'
#' @examples
#' # Setup model matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' colnames(design) <- paste0("ng", c(50, 100))
#' # Normalize data
#' yeast_norm <- yeast %>%
#'     psrn("identifier") %>%
#'     # Get mean-variance trends
#'     calculate_mean_sd_trends(design)
calculate_mean_sd_trends <- function(data, design_matrix){
  conditions <- design_matrix %>%
    colnames() %>%
    paste0(collapse = '|')
  data %>%
    dplyr::mutate(
      mean = rowMeans(dplyr::across(dplyr::matches(conditions)), na.rm = T),
      sd = apply(dplyr::across(dplyr::matches(conditions)), 1, sd, na.rm = T)
    )
}
