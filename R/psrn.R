utils::globalVariables(c("where", "value", "ref", "all_of"))
#' Normalize data to a pseudo-reference
#'
#' This function generates a pseudo-reference by taking the geometric mean of
#' each feature across all samples. Each feature in each sample is then divided
#' by the pseudo-reference. Then, the median ratio of all ratios is used as an
#' estimate to use for normalizing differences in loading concentration. All
#' features in each sample is then divided by their corresponding estimate.
#' All estimates are based on features without missing values.
#' For details see \insertCite{anders2010differential;textual}{baldur}.
#'
#' @param data data.frame
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#' @param load_info logical; should the load information be output?
#' @param log boolean variable indicating if the data should be log transformed
#'     after normalization
#' @param target target columns to normalize, supports
#'     \code{\link[tidyselect]{tidyselect-package}} syntax. By default, all numerical
#'     columns will be used in the normalization if not specified.
#'
#' @return data frame with normalized values if `load_info=FALSE`, if it is `TRUE`
#'    then it returns a list with two tibbles. One tibble containing the
#'    normalized data and one containing the loading info as well as the
#'    estimated normalization factors.
#' @export
#'
#' @importFrom dplyr .data
#' @importFrom dplyr enquo
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyr pivot_longer
#' @importFrom tibble enframe
#' @importFrom utils globalVariables
#'
#' @examples
#' yeast_psrn <- psrn(yeast, "identifier")
#' yeast_psrn_with_load <- psrn(yeast, "identifier", load_info = TRUE)
#' yeast_ng50_only <- psrn(yeast, "identifier", target = matches('ng50'))
#' @source \url{https://www.nature.com/articles/npre.2010.4282.1}
#' @references
#' \insertAllCited{}
psrn <- function(data,
                 id_col,
                 log = TRUE,
                 load_info = FALSE,
                 target = NULL) {
  check_id_col(id_col)
  target <- dplyr::enquo(target)
  target <- check_target(target)
  data_filtered <- data %>%
    tidyr::drop_na(!!target)
  loading_sizes <- calc_loading_size(data_filtered, target)
  pseudo_reference <- data_filtered %>%
    tidyr::pivot_longer(!!target) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      ref = prod(value^(1 / dplyr::n()))
    )
  scaling_factors <- data_filtered %>%
    tidyr::pivot_longer(!!target, names_to = "sample") %>%
    dplyr::left_join(pseudo_reference, by = id_col) %>%
    dplyr::mutate(
      value = value / ref
    ) %>%
    dplyr::select(-ref) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(rle_factor = stats::median(value)) %>%
    dplyr::left_join(loading_sizes, by = "sample")
  for (i in seq_len(nrow(scaling_factors))) {
    data[scaling_factors$sample[i]] <-
      data[scaling_factors$sample[i]] / scaling_factors$rle_factor[i]
  }
  if (log) {
    data <- data %>%
      dplyr::mutate(
        dplyr::across(!!target, log2)
      )
  }
  if (load_info) {
    return(
      list(
        data = data,
        scaling_factors = scaling_factors
      )
    )
  } else {
    return(data)
  }
}

calc_loading_size <- function(data, targets) {
  data %>%
    dplyr::select(!!targets) %>%
    colSums() %>%
    tibble::enframe(name = "sample", value = "load_size")
}

check_id_col <- function(id_col) {
  stopifnot(is.character(id_col))
  if (!id_col %in% names(data)) {
    stop(cat(id_col, '(id_col) is not a column in the dataset.'))
  }
}
