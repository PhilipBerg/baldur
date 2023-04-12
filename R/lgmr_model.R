#' Latent Gamma Regression Model
#'
#' @description Baldur has the option to use the Latent Gamma Regression Model
#'   (LGMR). This model attempts to estimate local Mean-Variance (M-V) trend for
#'   each peptide (or protein, PTM, etc.). To this end is describes the M-V
#'   trend as a regression with a common link function and a latent one. While
#'   there may not be a latent trend in real data, it allows for a mathematical
#'   description that can estimate the local trends of each peptide. The model
#'   describes the standard deviation (S) as a gamma random variable with mean
#'   dependency described by the sample mean (\eqn{\bar{y}}):
#'   \deqn{\boldsymbol{S}\sim\Gamma(\alpha,\beta),\quad\boldsymbol{\beta}=\frac{\alpha}{\boldsymbol{\mu}}}
#'   \deqn{\boldsymbol{\mu}=\exp(\gamma_0-\gamma_{\bar{y}}\,f(\boldsymbol{\bar{y}}))
#'   + \kappa
#'   \exp(\boldsymbol{\theta}(\gamma_{0L}-\gamma_{\bar{y}L}\,f(\boldsymbol{\bar{y}})),\quad
#'   f(\boldsymbol{x})=\frac{\boldsymbol{x}-\mu_{\bar{y}}}{\sigma_{\bar{y}}}}
#'   Here, \eqn{\exp(\gamma_0-\gamma_{\bar{y}})} is the common trend of the
#'   regression. Then, the mixing variable, \eqn{\boldsymbol{\theta}}, describes
#'   the proportion of the mixture that each peptide has, and \eqn{\kappa} is
#'   just some small constant such that when \eqn{\theta} is zero the latent
#'   term is small.
#'
#' @details Next we will describe the hierarchical prior setup for the
#' regression variables. A priori, Baldur assumes that each peptide is most
#' likely to completely have the mixture or to not have it at all. To this end,
#' the \eqn{\theta_i} for the i:th peptide has a beta distribution:
#' \deqn{\theta_i\sim\beta(0.5,0.5)} Further, it can be seen that Baldur assumes
#' that S always decreases (on average) with the sample mean. To this end, both
#' slopes (\eqn{\gamma_{\bar{y}},\gamma_{\bar{y}L}}) are defined on the real
#' positive line. Hence, we used Half-Normal (HN)
#' distributions to describe the slope parameters: \deqn{\gamma_{\bar{y}}\sim
#' HN(1)} \deqn{\gamma_{\bar{y}L}\sim HN(0.1)} For the intercepts, we assume a
#' standard normal prior for the common intercept. Then, we use a skew-normal to
#' model the latent intercept. The reason for this is two-fold, one,
#' \eqn{\kappa} will decrease the value of the latent term and, two, to push the
#' latent trend upwards. The latter reason is so that the latent intercept is
#' larger than the common and so that the priors prioritize a shift in intercept
#' over a increase in slope.
#' For the intercepts, Baldur uses a standard normal prior for the common intercept.
#' \deqn{\gamma_0\sim\mathcal{N}(0,1)}
#' While for the latent trend, it uses a skew-normal (SN) to push up the second
#' trend and to counteract the shrinkage of \eqn{\kappa}.
#' \deqn{\gamma_{0L}\sim\mathcal{SN}(2,15,35)}
#'
#'
#' @section Code: The `Stan` code for this model is given by:
#' ```{r, tidy = TRUE, comment = NA}
#' lgmr_model
#' ```
#'
#' @export
#' @name lgmr_model
#' @return A `stanmodel` that can be used in [fit_lgmr].
NULL
