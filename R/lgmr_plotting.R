utils::globalVariables(c("theta", "data", "x0", "x1", "y0", "y1", "dydx"))
#'Visualization of LGMR models
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'  Options to plot the LGMR model. `plot_lgmr_regression` will plot the data
#'  colored by the amount of latent trend they have as well as the two extreme
#'  regression cases when \eqn{\theta} is zero or one. `plot_regression_field`
#'  will plot the local regression trend for each data point as a vector field and
#'  color the vector based on the derivative at the mean of the peptide.
#'
#'@name lgmr_plotting
#'@param model An LGMR model object
#'
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 aes
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 theme_classic
#'@importFrom ggplot2 scale_color_viridis_c
#'@importFrom ggplot2 theme
#'@importFrom ggplot2 stat_function
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 geom_segment
#'@importFrom dplyr mutate
#'@importFrom rlang quo
#'@importFrom tidyr unnest
#'@importFrom purrr map
#'
#'
#'@return A ggplot object
#'@export
#'
#' @examples
#' #' # Define design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#'\donttest{
#' # Normalize data, calculate M-V trend, and fit LGMR model
#' yeast_lgmr <- yeast %>%
#'     # Remove missing values
#'     tidyr::drop_na() %>%
#'     # Normalize
#'     psrn("identifier") %>%
#'     # Add the mean-variance trends
#'     calculate_mean_sd_trends(design) %>%
#'     # Fit the model
#'     fit_lgmr("identifier")
#' # Print everything except thetas
#' print(yeast_lgmr, pars = c("coefficients", "auxiliary"))
#' # Extract the mean of the model parameters posterior
#' yeast_lgmr_pars <- coef(yeast_lgmr, pars = 'all', simplify = TRUE)
#' plot_lgmr_regression(yeast_lgmr)
#' plot_regression_field(yeast_lgmr)
#' }
plot_lgmr_regression <- function(model) {

  pars <- coef(model, TRUE, c('coef', 'theta'))
  reg_pars <- pars$coef

  mu_inputs <- mu_std_inputs(model$data)

  upper <- ~ mu_fun(rep(0, times = length(.x)),   reg_pars, .x, mu_inputs[1], mu_inputs[2])
  middl <- ~ mu_fun(rep(0.5, times = length(.x)), reg_pars, .x, mu_inputs[1], mu_inputs[2])
  lower <- ~ mu_fun(rep(1, times = length(.x)),   reg_pars, .x, mu_inputs[1], mu_inputs[2])

  model$data %>%
    dplyr::mutate(
      theta = pars$theta
    ) %>%
    ggplot2::ggplot(ggplot2::aes(mean, sd, color = theta)) +
    ggplot2::geom_point(size = .1) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_viridis_c(end = .9) +
    ggplot2::theme(
      legend.position = c(.75, .75)
    ) +
    ggplot2::stat_function(aes(color = 1),   fun = lower, linewidth = .9, n = 10000) +
    ggplot2::stat_function(aes(color = .5),  fun = middl, linewidth = .9, n = 10000) +
    ggplot2::stat_function(aes(color = 0),   fun = upper, linewidth = .9, n = 10000) +
    ggplot2::labs(
      x = expression(bold(bar(y))),
      y = expression(bold(s))
    )
}

#' @rdname lgmr_plotting
#'
#' @param n Number of interpolation points for each peptides regression
#' @param rng The proportional range of each peptides regression.
#'   E.g., a value of 10 will make each line span 1 % of the x-axis.
#'
#' @export
plot_regression_field <- function(model, n = 10, rng = 10) {
  rng <- diff(range(model$data$mean))/rng

  mu_inputs <- mu_std_inputs(model$data)

  pars <- coef(model, TRUE, c('coef', 'theta'))
  reg <- pars$coef

  deriv <- rlang::quo(
    (1/mu_inputs[2]) * (
      -reg['S_L']*.001*theta*exp(theta*(reg['I_L']-reg['S_L']*(mean - mu_inputs[1])/mu_inputs[2])) +
      -reg['S']*exp(reg['I']-reg['S']*(mean - mu_inputs[1])/mu_inputs[2])
    )
  )
  model$data %>%
    dplyr::mutate(
      theta = pars$theta,
      x0 = purrr::map(mean, ~seq(from = .x - rng, to = .x + rng, length.out = n)[1:(n - 1)]),
      x1 = purrr::map(mean, ~seq(from = .x - rng, to = .x + rng, length.out = n)[2:n])
    ) %>%
    tidyr::unnest(c(x0, x1)) %>%
    dplyr::mutate(
      y0 = mu_fun(theta, reg, x0, mu_inputs[1], mu_inputs[2]),
      y1 = mu_fun(theta, reg, x1, mu_inputs[1], mu_inputs[2]),
      dydx =  !!deriv
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x0, y0, xend = x1, yend = y1, color = dydx)) +
    ggplot2::geom_segment() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = expression(bold(bar(y))),
      y = expression(bold(s))
    ) +
    ggplot2::scale_color_viridis_c('Slope',
                          labels = ~ paste(ifelse(sign(.x) == -1, "\u2212", ""),
                                           sprintf("%2.2f", abs(.x))),
                          option = 'H', direction = -1
    )
}
