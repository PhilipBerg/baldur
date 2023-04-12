#' Baldur's weakly informative prior for the mean in conditions
#'
#' @description Here we will model the mean of the prior with a weakly
#'   informative (WI) prior. We will assume that, in essence, nothing is know
#'   about the mean. As such, for the WI prior, we use a normal prior on
#'   \eqn{\boldsymbol{\mu}_0} centered at zero and with a very large variance.
#'   \deqn{\boldsymbol{\mu}_0\sim\mathcal{N}(0,100)}
#'
#' @section Code: The `Stan` code for this model is given by:
#' ```{r, tidy = TRUE, comment = NA}
#' weakly_informative
#' ```
#'
#' @return A `stanmodel` that can be used in [infer_data_and_decision_model].
#' @export
#' @name weakly_informative
NULL
