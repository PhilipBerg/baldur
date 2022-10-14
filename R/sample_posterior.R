utils::globalVariables(c("alpha", "betau", "id", "tmp", "intu", "condi"))
#' Sample the Posterior of the Bayesian model
#'
#' @description Function to sample the posterior of the Bayesian decision model.
#'
#' @param data A `tibble` or `data.frame` with alpha,beta priors annotated
#' @param id_col_name A character of the id column
#' @param design_matrix A design matrix for the data (see example)
#' @param contrast_matrix A contrast matrix of the decisions to test such that the first column is the index of the column in the design matrix that should be compared to the second column
#' @param uncertainty_matrix A matrix of observation uncertainty
#' @param bayesian_model Which Bayesian model to use. Currently only one internal model allowed, this argument is for forward compatibility
#' @param clusters The number of parallel threads to run on.
#' @param robust If integration of the posterior should be done robust by trimming the tails
#' @param perc  If robust is `TRUE` how much of each tail to remove before integration. Should be between 0 (no trimming) and 0.5 (trimming all data).
#' @param mu_not The value of the null hypothesis for the difference in means
#' @param ... Additional arguments passed to \code{rstan::\link[rstan:sampling]{sampling}}. Note that verbose will always be forced to `FALSE`
#'
#' @return A `tibble` or `data.frame` annotated  with statistics of the posterior and p(error).
#' `err` is the probability of error, i.e., the two tail-density supporting the null-hypothesis, `lfc` is the estimated log$_2$-fold change, `sigma` is the common #' variance, and `lp` is the log-posterior.
#' Columns without suffix shows the mean estimate from the posterior, while the suffixes `_025`, `_50`, and `_975`, are the 2.5, 50.0, and 97.5, percentiles, respectively.
#' The suffixes `_eff` and `_rhat` are the diagnostic variables returned by `rstan` (please see the Stan manual for more details).
#' In general, a larger `_eff` indicates a better sampling efficiency, and `_rhat` compares the mixing within chains against between the chains and should be smaller than 1.05.
#' @export
#'
#' @importFrom rlang :=
#' @importFrom rlang dots_list
#' @importFrom rlang expr
#' @importFrom rlang sym
#' @importFrom rlang abort
#' @importFrom rlang current_env
#' @importFrom multidplyr new_cluster
#' @importFrom multidplyr cluster_library
#' @importFrom multidplyr cluster_copy
#' @importFrom multidplyr cluster_assign_partition
#' @importFrom multidplyr cluster_send
#' @importFrom multidplyr cluster_call
#' @importFrom multidplyr cluster_assign
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom dplyr transmute
#' @importFrom dplyr between
#' @importFrom dplyr filter
#' @importFrom dplyr tibble
#' @importFrom stringr word
#' @importFrom stats setNames
#' @importFrom stats ecdf
#' @importFrom stats quantile
#' @importFrom purrr map
#' @importFrom purrr pmap
#' @importFrom purrr map2_dfr
#' @importFrom purrr map_dbl
#' @importFrom purrr map2_chr
#' @importFrom purrr quietly
#' @importFrom tidyr unnest
#' @importFrom rstan extract
#' @importFrom rstan summary
#' @importFrom magrittr use_series
#'
#' @examples # (Please see the vignette for a tutorial)
#' # Setup model matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast_norm <- yeast %>%
#' # Remove missing data
#'   tidyr::drop_na() %>%
#'   # Normalize data
#'   psrn('identifier') %>%
#'   # Add mean-variance trends
#'   calculate_mean_sd_trends(design)
#' # Fit the gamma regression
#' gam <- fit_gamma_regression(yeast_norm, sd ~ mean)
#' # Estimate each data point's uncertainty
#' unc <- estimate_uncertainty(yeast_norm, 'identifier', design, gam)
#' yeast_norm <- yeast_norm %>%
#'    # Add hyper-priors for sigma
#'    estimate_gamma_priors(design, gam)
#' # Setup contrast matrix
#' contrast <- matrix(1:2, ncol = 2)
#' \donttest{
#' yeast_norm %>%
#'   sample_posterior(
#'     'identifier',
#'     design,
#'     contrast,
#'     unc,
#'     clusters = 1 # I highly recommend increasing the number of parallel workers/clusters
#'                  # this will greatly reduce the speed of running baldur
#'   )
#' }
sample_posterior <- function(data, id_col_name, design_matrix, contrast_matrix, uncertainty_matrix, bayesian_model = stanmodels$uncertainty_model_semi_informative, clusters = 1, robust = FALSE, perc = .05, mu_not = 0, ...){
  rstan_inputs <- rlang::dots_list(...)
  rstan_inputs$verbose <- F
  N <- sum(design_matrix)
  K <- ncol(design_matrix)
  C <- nrow(contrast_matrix)
  sigma_bar <- sd(data$mean)
  ori_data <- data
  pmap_columns <- rlang::expr(list(!!rlang::sym(id_col_name), alpha, beta))
  if(clusters != 1){
    cl <- multidplyr::new_cluster(clusters)
    suppressWarnings(
      invisible(
        purrr::quietly(multidplyr::cluster_library)(cl,
                                                    c(
                                                      "dplyr",
                                                      "tidyr",
                                                      "purrr",
                                                      "tibble",
                                                      "stringr",
                                                      "magrittr",
                                                      "rstan",
                                                      'StanHeaders',
                                                      'rlang'
                                                    )
        )
      )
    )
    multidplyr::cluster_copy(cl,
                             c(
                               "bayesian_testing",
                               "id_col_name",
                               "design_matrix",
                               "uncertainty_matrix",
                               "generate_stan_data_input",
                               "stan_summary",
                               "estimate_error",
                               "pmap_columns",
                               "ori_data",
                               "N",
                               "K",
                               "C",
                               "contrast_matrix",
                               'bayesian_model',
                               "robust",
                               "perc",
                               "mu_not",
                               "rstan_inputs",
                               "stan_nse_wrapper",
                               "sigma_bar"
                             )
    )
    model_check <- stringr::word(deparse(substitute(bayesian_model)), 2, sep = '\\$')
    if(!is.na(model_check)){
      multidplyr::cluster_copy(cl,
                               paste0(
                                 "rstantools_model_",
                                 stringr::word(deparse(substitute(bayesian_model)), 2, sep = '\\$')
                               )
      )
    }
    multidplyr::cluster_assign_partition(cl, id = ori_data$identifier, alpha = ori_data$alpha, beta = ori_data$beta)
    multidplyr::cluster_assign(cl,
                               condi = colnames(design_matrix),
                               condi_regex = paste0(colnames(design_matrix), collapse = '|')
    )
    multidplyr::cluster_send(cl,
                             dat <- purrr::pmap(list(stats::setNames(id, id), alpha, beta),
                                                ~ generate_stan_data_input(
                                                  ..1,
                                                  id_col_name,
                                                  design_matrix,
                                                  ori_data,
                                                  uncertainty_matrix,
                                                  contrast_matrix,
                                                  N, K, C,
                                                  ..2, ..3,
                                                  condi,
                                                  condi_regex,
                                                  sigma_bar
                                                )
                             )
    )
    results <- multidplyr::cluster_call(cl,
                                        purrr::map(dat,
                                                   stan_nse_wrapper,
                                                   bayesian_model,
                                                   rstan_inputs
                                        ) %>%
                                          purrr::map2_dfr(dat,
                                                          stan_summary,
                                                          condi,
                                                          C,
                                                          robust,
                                                          perc,
                                                          mu_not,
                                                          .id = id_col_name
                                          )
    ) %>%
      dplyr::bind_rows()
    rm(cl)
    gc(FALSE)
    return(results)
  } else {
    data %>%
      dplyr::transmute(
        !!id_col_name := .data[[id_col_name]],
        tmp = purrr::pmap(!!pmap_columns,
                          bayesian_testing,
                          ori_data, id_col_name, design_matrix, contrast_matrix, bayesian_model, N, K, C, uncertainty_matrix, robust, perc, mu_not, rstan_inputs, sigma_bar
        )
      ) %>%
      tidyr::unnest(tmp)
  }
}

stan_nse_wrapper <- function(data, model, ...){
  rlang::eval_tidy(
    rlang::call2(
      rstan::sampling,
      object = model,
      data = data,
      !!!rlang::dots_splice(...)
    )
  )
}

generate_stan_data_input <- function(id, id_col_name, design_matrix, data, uncertainty, comparison, N, K, C, alpha, beta, condi, condi_regex, sigma_bar){
  row <- data[data[[id_col_name]] == id, stringr::str_detect(names(data), condi_regex)]
  ybar <- purrr::map_dbl(purrr::set_names(condi, condi),
                         ~as.numeric(row[stringr::str_which(colnames(row), .x)]) %>%
                           mean()
  )[condi]
  if(is.null(uncertainty)){
    u = rep(1,N)
  }else{
    u = uncertainty[id,]
  }
  list(
    N = N,
    K = K,
    C = C,
    x = design_matrix,
    y = as.numeric(row),
    c = comparison,
    alpha = alpha,
    beta_gamma = beta,
    mu_not = ybar,
    u = u,
    sigma_mu_not = sigma_bar
  )
}

estimate_error <- function(posterior, robust, perc, mu_not) {
  if (robust) {
    if (perc >= .5) {
      rlang::abort(
        c(
          'perc must be less than 0.5 when robust = TRUE or all posetrior data will be filtered.',
          'Please set perc to a value less than 0.5, we do not recommend larger than 0.25 and a large sample from the posterior would be needed.'
        ),
        class = 'posterior-filtering',
        call = rlang::current_env()
      )
    }
    q <- stats::quantile(posterior, c(perc, 1-perc))
    posterior <- posterior[ dplyr::between(posterior, q[1], q[2]) ]
  }
  p <- pnorm(mu_not, mean(posterior), sd(posterior))
  if (p < .5) {
    2 * p
  }
  else {
    2 * (1 - p)
  }
}

stan_summary <- function(fit, dat, condi, C, robust, perc, mu_not){
  lfc_pars <- paste0('y_diff[', seq_len(C), ']')
  err <- rstan::extract(fit, pars = lfc_pars) %>%
    purrr::map_dbl(estimate_error, robust, perc, mu_not)
  summ <- rstan::summary(fit, pars = c(lfc_pars, 'sigma', 'lp__')) %>%
    magrittr::use_series(summary)
  summ <- summ[
    rownames(summ) %in% c(lfc_pars, 'sigma', 'lp__'),
    colnames(summ) %in% c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat', '50%')
  ]
  summ <- summ[1:nrow(summ), c(1:2, 4:6, 3)]
  dplyr::tibble(
    comparison = purrr::map2_chr(condi[dat$c[,1]], condi[dat$c[,2]], ~stringr::str_flatten(c(.x, .y), collapse = ' vs ')),
    err = err,
    lfc = summ[rownames(summ) %in% lfc_pars, 1],
    lfc_025 = summ[rownames(summ) %in% lfc_pars, 2],
    lfc_50 = summ[rownames(summ) %in% lfc_pars, 6],
    lfc_975 = summ[rownames(summ) %in% lfc_pars, 3],
    lfc_eff = summ[rownames(summ) %in% lfc_pars, 4],
    lfc_rhat = summ[rownames(summ) %in% lfc_pars, 5],
    sigma = summ[rownames(summ) == 'sigma', 1],
    sigma_025 = summ[rownames(summ) == 'sigma', 2],
    sigma_50 = summ[rownames(summ) == 'sigma', 6],
    sigma_975 = summ[rownames(summ) == 'sigma', 3],
    sigma_eff = summ[rownames(summ) == 'sigma', 4],
    sigma_rhat = summ[rownames(summ) == 'sigma', 5],
    lp = summ[rownames(summ) == 'lp__', 1],
    lp_025 = summ[rownames(summ) == 'lp__', 2],
    lp_50 = summ[rownames(summ) == 'lp__', 6],
    lp_975 = summ[rownames(summ) == 'lp__', 3],
    lp_eff = summ[rownames(summ) == 'lp__', 4],
    lp_rhat = summ[rownames(summ) == 'lp__', 5]
  )
}

bayesian_testing <- function(id, alpha, beta, data, id_col_name, design_matrix, comparison, model, N, K, C, uncertainty = NULL, robust, perc, mu_not, rstan_args, sigma_bar){
  condi <- colnames(design_matrix)
  condi_regex <- paste0(condi, collapse = '|')
  dat <- generate_stan_data_input(id, id_col_name, design_matrix, data, uncertainty, comparison, N, K, C, alpha, beta, condi, condi_regex, sigma_bar)
  purrr::quietly(stan_nse_wrapper)(dat, model, rstan_args)$result %>%
    stan_summary(dat, condi, C, robust, perc, mu_not)
}
