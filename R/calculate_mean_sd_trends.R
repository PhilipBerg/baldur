#' Calculate the Mean-Variance trend
#'
#' @description Calculates the mean and standard deviation of each row (peptide)
#'   and adds them as new columns. Assumes that the condition names are the
#'   names in the design matrix.
#' @param data A `tibble` or `data.frame` to annotate with mean and sd
#' @param design_matrix A design matrix for the data (see example).
#'
#' @return A `tibble` or `data.frame` with the mean and sd vectors
#' @export
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr matches
#'
#' @examples
#'
#' # Setup model matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast %>%
#'     # Normalize data
#'     psrn("identifier") %>%
#'     # Get mean-variance trends
#'     calculate_mean_sd_trends(design)
calculate_mean_sd_trends <- function(data, design_matrix){
  check_design(design_matrix, data)

  conditions <- design_matrix %>%
    colnames() %>%
    paste0(collapse = '|')

  data %>%
    dplyr::mutate(
      mean = rowMeans(dplyr::across(dplyr::matches(conditions)), na.rm = T),
      sd = apply(dplyr::across(dplyr::matches(conditions)), 1, sd, na.rm = T)
    )
}

check_design <- function(design_matrix, data, cenv = rlang::caller_call()) {
  rlang::is_missing(design_matrix)
  if (is.null(colnames(design_matrix))) {
    rlang::abort(
      c(
        'Design matrix has no column names',
        "Please see vignette('baldur_yeast_tutorial') or vignette('baldur_ups_tutorial') on how to setup the design matrix"
      ),
      call = cenv
    )
  }

  condi_counts <- colSums(design_matrix)
  col_counts <- map_dbl(names(condi_counts),
                        ~ sum(stringr::str_count(colnames(data), .x))
  )

  if (any(condi_counts != col_counts)) {
    rlang::abort(
      c(
        'Number of samples in design matrix does not match number of samples in data.',
        paste0(
          'Design matrix has: ', stringr::str_flatten(condi_counts, ', '),
          ' samples in conditions ', stringr::str_flatten(names(condi_counts), ', '), ', respectively'
        ),
        paste0(
          'Samples found in data are: ', stringr::str_flatten(col_counts, ', '),
          ' samples in conditions ', stringr::str_flatten(names(condi_counts), ', '), ', respectively'
        ),
        "They need to match, see vignette('baldur_yeast_tutorial') or vignette('baldur_ups_tutorial') on how to setup the design matrix"
      ),
      call = cenv
    )
  }
}
