#' @importFrom Rdpack reprompt
utils::globalVariables(c("sd", "model"))
#' Function for plotting the mean-variance gamma regressions
#'
#' Generates a scatter plot with the gamma regressions of the mean-variance
#' trend
#'
#' @param data The data to use for producing the plots.
#' @param design A design matrix as produced by
#'   \code{\link[stats]{model.matrix}}.
#'
#' @return a plot with the mean-variance trend before partitioning on the left
#'   side, and the right side after.
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
#'     psrn("identifier") %>%
#'     plot_gamma_regression(design)
#'
plot_gamma_regression <- function(data, design) {
  data <- data %>%
    calculate_mean_sd_trends(design)
  data %>%
    plot_gamma()
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
    plot_mean_sd_trend()
}
