#' Title
#'
#' @name lgmr_plotting
#' @param model
#'
#' @return
#' @export
#'
#' @examples
plot_lgrm_regression <- function(model) {

  pars <- coef(model, TRUE, c('coef', 'theta'))
  reg_pars <- pars$coef

  mu_inputs <- mu_std_inputs(model$data)

  upper <- ~ mu_fun(0,   reg_pars, .x, mu_inputs[1], mu_inputs[2])
  middl <- ~ mu_fun(0.5, reg_pars, .x, mu_inputs[1], mu_inputs[2])
  lower <- ~ mu_fun(1,   reg_pars, .x, mu_inputs[1], mu_inputs[2])

  plt <- model$data %>%
    mutate(
      theta = pars$theta
    ) %>%
    ggplot(aes(mean, sd, color = theta)) +
    geom_point(size = .1) +
    theme_classic() +
    scale_color_viridis_c(end = .9) +
    theme(
      legend.position = c(.75, .75)
    ) +
    stat_function(aes(color = 1),   fun = lower, linewidth = .9, n = 10000) +
    stat_function(aes(color = .5),  fun = middl, linewidth = .9, n = 10000) +
    stat_function(aes(color = 0),   fun = upper, linewidth = .9, n = 10000) +
    labs(
      x = expression(bold(bar(y))),
      y = expression(bold(s))
    )
}

#' Vector Field Representation Of The LGMR Model
#'
#' @rdname lgmr_plotting
#'
#' @param model
#' @param n
#' @param rng
#'
#' @return
#' @export
#'
#' @examples
plot_regression_field <- function(model, n = 10, rng = NULL) {
  if(is.null(rng)) {
    rng <- diff(range(data$mean))/10
  } else if (is.numeric(rng)) {
    rng <- diff(range(data$mean))/rng
  }

  mu_inputs <- mu_std_inputs(model$data)

  pars <- coef(model, TRUE, c('coef', 'theta'))
  reg <- pars$coef
  theta <- pars$theta

  deriv <- quo(
    -reg['S_L']*.001*theta*exp(theta*(reg['I_L']-reg['S_L']*(mean - mu_inputs[1])/mu_inputs[2])) +
      -reg['S']*exp(reg['I']-reg['S']*(mean - mu_inputs[1])/mu_inputs[2])
  )
  model$data %>%
    mutate(
      x0 = map(mean, ~seq(from = .x - rng, to = .x + rng, length.out = n)[1:(n - 1)]),
      x1 = map(mean, ~seq(from = .x - rng, to = .x + rng, length.out = n)[2:n])
    ) %>%
    unnest(c(x0, x1)) %>%
    mutate(
      y0 = mu_fun(theta, reg, x0, m, s),
      y1 = mu_fun(theta, reg, x1, m, s),
      dydx =  !!deriv
    ) %>%
    ggplot(aes(x0, y0, xend = x1, yend = y1, color = dydx)) +
    geom_segment() +
    theme_classic() +
    labs(
      x = expression(bold(bar(y))),
      y = expression(bold(s))
    ) +
    scale_color_viridis_c('Slope',
                          labels = ~ paste(ifelse(sign(.x) == -1, "\u2212", ""),
                                           sprintf("%2.2f", abs(.x))),
                          option = 'H', direction = -1
    )
}
