#' @importFrom Rdpack reprompt
utils::globalVariables(c(".", "sd", "model"))
#' Function for plotting the mean-variance gamma regressions
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends with and without partitioning
#'
#' @param data The data to use for producing the plots.
#' @param design A design matrix as produced by \code{\link[stats]{model.matrix}}.
#' @param ... Additional arguments to \code{\link[baldur]{trend_partitioning}}.
#'
#' @return a plot with the mean-variance trend before partitioning on the left side, and the right side after.
#' @export
#'
#' @importFrom cowplot plot_grid
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_label
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 stat_function
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 theme
#' @importFrom purrr set_names
#' @importFrom rlang call2
#' @importFrom rlang eval_tidy
#' @importFrom rlang dots_list
#' @importFrom stats Gamma
#' @importFrom stats glm
#' @importFrom stats predict.glm
#' @importFrom tidyr drop_na
#' @importFrom viridisLite turbo
#'
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast_norm <- yeast %>%
#'     # Remove missing data
#'     # Note that this could be replaced with imputation
#'     tidyr::drop_na() %>%
#'     # Normalize
#'     psrn("identifier")
#'
#' # Generate the plots
#' \donttest{
#' plot_gamma_regression(yeast_norm, design, verbose = FALSE)
#' }
plot_gamma_regression <- function(data, design, ...) {
  if(!'sd' %in% names(data)){
    data <- data %>%
      calculate_mean_sd_trends(design)
  }
  base_plot <- data %>%
    plot_gamma()
  part_plot <- data %>%
    plot_gamma_partition(design, ...)
  plots <- cowplot::plot_grid(base_plot, part_plot)
  title <- cowplot::ggdraw() +
    cowplot::draw_label("Mean-Variance trends", fontface = "bold")
  cowplot::plot_grid(
    title,
    plots,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
}

plot_mean_sd_trend <- function(data) {
  data %>%
    ggplot2::ggplot(ggplot2::aes(mean, sd)) +
    ggplot2::geom_point(size = 1 / 10) +
    ggplot2::geom_smooth(
      method = stats::glm,
      formula = y ~ x,
      method.args = list(family = stats::Gamma(log)),
      fullrange = TRUE,
      se = F,
      color = 'grey'
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = expression(bold(bar(y))), y = expression(bold(s))
    )
}

#' Function for plotting the gamma regression for the mean-variance trend
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends without partitioning
#'
#' @param data The data to use for producing the plots.
#'
#' @return a plot with the estimated mean-variance trend
#' @export
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast_norm <- psrn(yeast, "identifier")
#'
#' # Generate the plot
#' yeast_norm %>%
#'   calculate_mean_sd_trends(design) %>%
#'   plot_gamma()
plot_gamma <- function(data) {
  data %>%
    tidyr::drop_na(sd) %>%
    plot_mean_sd_trend() +
    ggplot2::ggtitle("Before Partitioning")
}

#' Function for plotting the gamma regression after partitioning
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trends after partitioning.
#'
#' @param data The data to use for producing the plots.
#' @param design A design matrix as produced by \code{stats::\link[stats:model.matrix]{model.matrix}}.
#' @param ... Additional arguments to \code{baldur::\link[baldur:trend_partitioning]{trend_partitioning}}.
#'
#' @return a plot with the mean-variance trend after partitioning
#' @export
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast_norm <- yeast %>%
#'     # Remove missing data
#'     # Note that this could be replaced with imputation
#'     tidyr::drop_na() %>%
#'     # Normalize
#'     psrn("identifier")
#'
#' # Generate the plot
#' \donttest{
#' yeast_norm %>%
#'   calculate_mean_sd_trends(design) %>%
#'   plot_gamma_partition(design, verbose = FALSE)
#' }
plot_gamma_partition <- function(data, design, ...) {
  if(!'c' %in% names(data)){
    data <- data %>%
      trend_partitioning(design, ...)
  }
  trend_colors <- purrr::set_names(viridisLite::turbo(2, end = .75), c('Lower', 'Upper'))
  gam_reg <- rlang::eval_tidy(
    rlang::call2(fit_gamma_regression, data = data, !!!rlang::dots_list(...))
  )
  data %>%
    dplyr::mutate(
      c = stringr::str_replace_all(c, purrr::set_names(c('Lower', 'Upper'), c('L', 'U')))
    ) %>%
    tidyr::drop_na(sd) %>%
    ggplot2::ggplot(ggplot2::aes(mean, sd, color = c)) +
    ggplot2::geom_point(size = 1 / 10) +
    ggplot2::theme_classic() +
    ggplot2::stat_function(
      fun = ~stats::predict.glm(gam_reg, newdata = data.frame(mean = .x, c = 'L'), type = 'response'), color = 'blue'
    ) +
    ggplot2::stat_function(
      fun = ~stats::predict.glm(gam_reg, newdata = data.frame(mean = .x, c = 'U'), type = 'response'), color = 'red'
    ) +
    ggplot2::ggtitle("After Partitioning") +
    ggplot2::scale_color_manual('Trend', values = trend_colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2))) +
    ggplot2::labs(
      x = expression(bold(bar(y))), y = expression(bold(s))
    ) +
    ggplot2::theme(
      legend.position = c(.8, .9)
    )
}
