#' Baldur's empirical Bayes Prior For The Mean In Conditions
#'
#' @description
#' Here we assume that the sample mean of each condition is an estimator for the center of the mean prior.
#' In addition, it assumes that the confidence in the prior is proportional to the variance of the peptide.
#' \deqn{\boldsymbol{\mu}_0\sim\mathcal{N}(\boldsymbol{\bar{y}},\sigma\boldsymbol{n}_R)}
#' \deqn{\boldsymbol{n}_R=[\frac{1}{\sqrt{n_1}},\frac{1}{\sqrt{n_2}},\dots,\frac{1}{\sqrt{n_C}}]}
#'
#' @section Code: The `Stan` code for this model is given by:
#' ```{r, tidy = TRUE, comment = NA}
#' empirical_bayes
#' ```
#'
#' @return A `stanmodel` that can be used in [infer_data_and_decision_model].
#' @export
#' @name empirical_bayes
NULL
