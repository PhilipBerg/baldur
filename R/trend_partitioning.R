#' Mean-Variance Trend Partitioning
#'
#' @description Partitions the data into two mixtures.
#' @param data A `tibble` or `data.frame` to partition
#' @param design_matrix A design matrix for the data (see example)
#' @param formula Formula for the Gamma regression
#' @param eps Size of the integration window
#' @param n Number of integration windows
#' @param verbose If the number of points moved should be output
#'
#' @return A `tibble` or `data.frame` the partitioning vector `c`
#' @export
#'
#' @examples 'lorem'
trend_partitioning <- function(data, design_matrix, formula = sd ~ mean + c, eps = .1, n = 1000, verbose = T){
  cur_data <- data %>%
    prep_data_for_clustering(design_matrix, eps = eps, n = n) %>%
    run_procedure(formula, eps = eps, n = n)
  while (cur_data$i > 1) {
    if(verbose){
      print(
        paste('Moved', cur_data$i, 'points')
      )
    }
    cur_data <- cur_data$data %>%
      run_procedure(formula, eps = eps, n = n)
  }
  cur_data$data %>%
    dplyr::select(-c(betal, betau, intl, intu, alpha))
}

prep_data_for_clustering <- function(data, design_matrix, eps = .1, n = 1000){
  data_ms <- data %>%
    calculate_mean_sd_trends(design_matrix)
  gam_reg <-  glm(sd ~ mean, Gamma(log), data_ms)
  data_ms %>%
    dplyr::mutate(
      res = residuals(gam_reg),
      c = dplyr::if_else(res<0, 'L', 'U')
    ) %>%
    dplyr::select(-res)
}

clust_itt_norm <- function(data, reg){
  tmp <- data$c
  data <- data %>%
    dplyr::mutate(
      c = dplyr::if_else(intl<intu, 'U', 'L')
    )
  i <- sum(tmp != data$c)
  return(
    list(
      data = data,
      i = i
    )
  )
}

resetimate_gamma_pars <- function(data, formula){
  gm <- glm(formula, Gamma(log), data)
  data %>%
    dplyr::mutate(
      alpha = (1/summary(gm)$dispersion),
      betal = estimate_beta(gm, mean, 'L', alpha),
      betau = estimate_beta(gm, mean, 'U', alpha)
    )
}

add_integrals <- function(data, eps = .1, n = 1000){
  data %>%
    dplyr::mutate(
      intu = purrr::pmap_dbl(list(sd, alpha, betau), num_int_trapz, eps, n),
      intl = purrr::pmap_dbl(list(sd, alpha, betal), num_int_trapz, eps, n)
    )
}

run_procedure <- function(data, formula, eps = .1, n = 1000){
  data %>%
    resetimate_gamma_pars(formula) %>%
    add_integrals(eps = eps, n = n) %>%
    clust_itt_norm()
}

num_int_trapz <- function(sd, alpha, beta, eps, n){
  rng <- c(-eps, eps) + sd
  x <- seq(rng[1], rng[2], length.out = n)
  dx <- diff(x)
  y <- dgamma(x, alpha, beta)
  sum(c(y[1]/2, y[2:(length(y)-2)], y[length(y)-1]/2)*dx)
}

estimate_beta <- function(model, mean, c, alpha){
  alpha / predict.glm(model, newdata = data.frame(mean = mean, c = c), type = 'response')
}
