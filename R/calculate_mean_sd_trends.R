#' Calculate the Mean-Variance trend
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#'  Calculates the mean and standard deviation of each row (peptide)
#'  and adds them as new columns. Assumes that the condition names are the
#'  names in the design matrix.
#'
#' @param data A `tibble` or `data.frame` to annotate with mean and sd
#' @param design_matrix A design matrix for the data (see example).
#' @param auxiliary_columns Names of columns in the design matrix that does not
#' have corresponding data in the data set. For example, this can be
#' co-founding variables such as bashes.
#'
#' @return A `tibble` or `data.frame` with the mean and sd vectors
#' @export
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr matches
#' @importFrom utils combn
#' @importFrom utils stack
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
calculate_mean_sd_trends <- function(data, design_matrix, auxiliary_columns = c()){
  design_matrix <- check_design_aux(design_matrix, auxiliary_columns)
  check_design(design_matrix, data)

  conditions <- get_conditions(design_matrix)

  data %>%
    dplyr::mutate(
      mean = rowMeans(dplyr::across(dplyr::matches(conditions)), na.rm = T),
      sd = apply(dplyr::across(dplyr::matches(conditions)), 1, sd, na.rm = T)
    )
}

check_design <- function(design_matrix, data, cenv = rlang::caller_call()) {
  inp <- rlang::enquo(design_matrix)
  rlang::is_missing(design_matrix)
  check_matrix(inp, "design matrix", cenv)
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
  condi <- names(condi_counts)
  cols <- data %>%
    dplyr::select(matches(stringr::str_flatten(condi, '|'))) %>%
    colnames()

  check_dups <- apply(utils::combn(condi, 2), 2, check_column_directionality, simplify = T)
  check_dups <- check_dups[!is.na(check_dups)]

  tbl <- setNames(paste0("^", condi), condi) %>%
    purrr::map(grep, cols, value = TRUE) %>%
    utils::stack() %>%
    table()
  dup_rows <- rowSums(tbl) != 1

  if(length(check_dups) != 0) {
    for (i in seq_along(check_dups)) {
      tbl <- correct_dups(check_dups[[1]], tbl, dup_rows)
    }
  }

  col_counts <- colSums(tbl)
  col_counts <- col_counts[names(condi_counts)]
  cols <- which(condi_counts != col_counts)

  if (length(cols) != 0) {
    rlang::abort(
      c(
        'Number of samples in design matrix does not match number of samples in data.',
        paste0(
          'Design matrix has: ', stringr::str_flatten(condi_counts[cols], ', '),
          ' samples in conditions ', stringr::str_flatten(names(condi_counts)[cols], ', '), ', respectively'
        ),
        paste0(
          'Samples found in data are: ', stringr::str_flatten(col_counts[cols], ', '),
          ' samples in conditions ', stringr::str_flatten(names(condi_counts)[cols], ', '), ', respectively'
        ),
        "They need to match, see vignette('baldur_yeast_tutorial') or vignette('baldur_ups_tutorial') on how to setup the design matrix."
      ),
      call = cenv
    )
  }
}

check_design_aux <- function(design_matrix, auxiliary_columns) {
  if (!is.null(auxiliary_columns)) {
    if (is.character(auxiliary_columns))
      auxiliary_columns <- which(colnames(design_matrix) %in% auxiliary_columns)
    design_matrix <- design_matrix[, -auxiliary_columns, drop = F]
  }
  return(design_matrix)
}

check_matrix <- function(input, type, cenv) {
  inp <- eval_tidy(input)
  inp_name <- rlang::as_name(input)
  if (!is.matrix(inp)) {
    rlang::abort(
      c(
        paste0(
          "Your ", type, ", ", inp_name, ", needs to be a matrix not a ",
          stringr::str_flatten(class(inp), ", "),
          "."
        )
      ),
      call = cenv
    )
  }
}

get_conditions <- function(design_matrix) {
  design_matrix %>%
    colnames() %>%
    paste0('^', ., collapse = '|')
}

check_column_directionality <- function(nams) {
  nams_reg <- paste0("^", nams)
  check <- c(
    grepl(nams_reg[1], nams[2]),
    grepl(nams_reg[2], nams[1])
  )
  if (any(check != 0)) {
    c(nams[check], nams[!check])
  } else{
    NA
  }
}

correct_dups <- function(dups, tab, rws) {
  tab[rowSums(tab[,dups]) == 2, dups[1]] <- 0
  return(tab)
}
