#' @importFrom rlang quo_is_null
#' @importFrom rlang expr
check_target <- function(target) {
  if (rlang::quo_is_null(target)) {
    rlang::expr(where(is.numeric))
  } else {
    target
  }
}
