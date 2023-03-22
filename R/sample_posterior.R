utils::globalVariables(c("alpha", "betau", "id", "tmp", "intu", "condi"))
#' Sample the Posterior of the data and decision model and generate point
#' estimates
#'
#' @description Function to sample the posterior of the Bayesian data and
#'   decision model. It first produces the needed inputs for Stan's [sampling()]
#'   for each peptide (or protein, PTM, etc.). It then runs the sampling for the
#'   data and decision model. From the posterior, it then collects estimates and
#'   sampling statistics from the posterior of data model and integrates the
#'   decision distribution D. It then returns a [tibble()] with all the
#'   information for each peptide's posterior (see **Value** below). There are
#'   major time gains to be made by running this procedure in parallel.
#'   `infer_data_and_decision_model` has an efficient wrapper around
#'   `multipldyr`. This will let you to evenly distribute all peptides evenly to
#'   each worker. E.g., two workers will each run half of the peptides in
#'   parallel.
#'
#' @details Model description here
#'
#' @param data A `tibble` or `data.frame` with alpha,beta priors annotated
#' @param id_col A character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#' @param design_matrix A design matrix for the data (see example)
#' @param contrast_matrix A contrast matrix of the decisions to test. Columns
#'   should sum to `0` and only mean comparisons are allowed. That is, the
#'   absolute value of the positive and negative values in each column has to
#'   sum to `2`. E.g., a column can be `[`0.5, 0.5, -1`]`\eqn{^T} but not `[`1,
#'   1, -1`]` or `[`0.5, 0.5, -2`]`. That is, `sum(abs(x))=2` where `x` is a
#'   column in the contrast matrix.
#' @param uncertainty_matrix A matrix of observation specific uncertainties
#' @param stan_model Which Bayesian model to use. Defaults to [empirical_bayes]
#'   but also allows [weakly_informative], or an user supplied function see [].
#' @param clusters The number of parallel threads/workers to run on.
#' @param h_not The value of the null hypothesis for the difference in means
#' @param ... Additional arguments passed to
#'   \code{\link[rstan:sampling]{rstan::sampling}}. Note that verbose will
#'   always be forced to `FALSE` to avoid console flooding.
#'
#' @return A [tibble()] or [data.frame()] annotated  with statistics of the
#'   posterior and p(error). `err` is the probability of error, i.e., the two
#'   tail-density supporting the null-hypothesis, `lfc` is the estimated
#'   \eqn{\log_2}-fold change, `sigma` is the common variance, and `lp` is the
#'   log-posterior. Columns without suffix shows the mean estimate from the
#'   posterior, while the suffixes `_025`, `_50`, and `_975`, are the 2.5, 50.0,
#'   and 97.5, percentiles, respectively. The suffixes `_eff` and `_rhat` are
#'   the diagnostic variables returned by `Stan`. In general, a larger `_eff`
#'   indicates effective sample size and `_rhat` indicates the mixing
#'   within chains and between the chains and should be smaller than 1.05
#'   (please see the Stan manual for more details).
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
#' @importFrom stats pnorm
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
#'
#' # Fit the gamma regression
#' gam <- fit_gamma_regression(yeast_norm, sd ~ mean)
#'
#' # Estimate each data point's uncertainty
#' unc <- estimate_uncertainty(gam, yeast_norm, 'identifier', design)
#'
#' yeast_norm <- gam %>%
#'    # Add hyper-priors for sigma
#'    estimate_gamma_hyperparameters(yeast_norm)
#' # Setup contrast matrix
#' contrast <- matrix(c(-1, 1)
#' \donttest{
#' yeast_norm %>%
#'   head() %>% # Just running a few for the example
#'   infer_data_and_decision_model(
#'     'identifier',
#'     design,
#'     contrast,
#'     unc,
#'     clusters = 1 # I highly recommend increasing the number of parallel workers/clusters
#'                  # this will greatly reduce the speed of running baldur
#'   )
#' }
infer_data_and_decision_model <- function(data, id_col, design_matrix, contrast_matrix, uncertainty_matrix, stan_model = empirical_bayes, clusters = 1, h_not = 0, ...){
  check_contrast_and_design(contrast_matrix, design_matrix)
  rstan_inputs <- rlang::dots_list(...)
  rstan_inputs$verbose <- F
  N <- sum(design_matrix)
  K <- ncol(design_matrix)
  C <- ncol(contrast_matrix)
  sigma_bar <- sd(data$mean)
  ori_data <- data
  pmap_columns <- rlang::expr(list(!!rlang::sym(id_col), alpha, beta))
  if (clusters > nrow(data)) {
    warning("You have more workers than rows in your data, setting workers to the number of rows in the data.")
    clusters <- nrow(data)
  }
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
                               "id_col",
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
                               'stan_model',
                               "h_not",
                               "rstan_inputs",
                               "stan_nse_wrapper",
                               "sigma_bar"
                             )
    )

    multidplyr::cluster_assign_partition(cl,
                                         id = ori_data$identifier,
                                         alpha = ori_data$alpha,
                                         beta = ori_data$beta
    )
    multidplyr::cluster_assign(cl,
                               condi = colnames(design_matrix),
                               condi_regex = paste0(colnames(design_matrix), collapse = '|'
                               )
    )
    multidplyr::cluster_send(cl,
                             dat <- purrr::pmap(list(stats::setNames(id, id),
                                                     alpha, beta),
                                                ~ generate_stan_data_input(
                                                  ..1,
                                                  id_col,
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
                                                   stan_model,
                                                   rstan_inputs
                                        ) %>%
                                          purrr::map2_dfr(dat,
                                                          stan_summary,
                                                          condi,
                                                          contrast_matrix,
                                                          h_not,
                                                          .id = id_col
                                          )
    ) %>%
      dplyr::bind_rows()
    rm(cl)
    gc(FALSE)
    return(results)
  } else {
    data %>%
      dplyr::transmute(
        !!id_col := .data[[id_col]],
        tmp = purrr::pmap(!!pmap_columns,
                          bayesian_testing,
                          ori_data, id_col, design_matrix, contrast_matrix,
                          stan_model, N, K, C, uncertainty_matrix, h_not,
                          rstan_inputs, sigma_bar
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

generate_stan_data_input <- function(id, id_col, design_matrix, data, uncertainty, comparison, N, K, C, alpha, beta, condi, condi_regex, sigma_bar){
  row <- data[data[[id_col]] == id, stringr::str_detect(names(data), condi_regex)]
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
    beta = beta,
    mu_not = ybar,
    u = u
  )
}

estimate_error <- function(posterior, h_not) {
  s  <- -abs(mean(posterior) - h_not)/sd(posterior)
  2*pnorm(s)
}

stan_summary <- function(fit, dat, condi, contrast,  h_not){
  lfc_pars <- paste0('y_diff[', seq_len(ncol(contrast)), ']')
  mu_pars  <- paste0("mu[", seq_along(condi), "]")
  err <- rstan::extract(fit, pars = lfc_pars) %>%
    purrr::map_dbl(estimate_error, h_not)
  summ <- rstan::summary(fit, pars = c(lfc_pars, mu_pars, 'sigma', 'lp__')) %>%
    magrittr::use_series(summary)
  summ <- summ[
    rownames(summ) %in% c(lfc_pars, mu_pars, 'sigma', 'lp__'),
    colnames(summ) %in% c('mean', '2.5%', '50%', '97.5%', 'n_eff', 'Rhat')
  ]
  summ <- summ[, c(1:2, 4:6, 3)]
  comps <- contrast_to_comps(contrast, condi)
  dplyr::tibble(
    comparison = comps,
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

bayesian_testing <- function(id, alpha, beta, data, id_col, design_matrix, comparison, model, N, K, C, uncertainty = NULL, h_not, rstan_args, sigma_bar){
  condi <- colnames(design_matrix)
  condi_regex <- paste0(condi, collapse = '|')
  dat <- generate_stan_data_input(id, id_col, design_matrix, data, uncertainty, comparison, N, K, C, alpha, beta, condi, condi_regex, sigma_bar)
  stan_output <- purrr::quietly(stan_nse_wrapper)(dat, model, rstan_args)
  stan_output$result %>%
    stan_summary(dat, condi, dat$c, h_not) %>%
    dplyr::mutate(
      warnings = list(stan_output$warnings)
    )
}

contrast_to_comps <- function(contrast, conditions){
  positives <- apply(contrast, 2, \(x) which(x>0)) %>%
    purrr::map_chr(~ stringr::str_flatten(conditions[.x], ' and '))
  negatives <- apply(contrast, 2, \(x) which(x<0)) %>%
    purrr::map_chr(~ stringr::str_flatten(conditions[.x], ' and '))
  stringr::str_c(positives, negatives, sep = ' vs ')
}

check_contrast_and_design <- function(contrast_matrix, design_matrix) {
  error <- FALSE
  cs  <- colSums(contrast_matrix)
  acs <- colSums(abs(contrast_matrix))
  cr <- nrow(contrast_matrix)
  dn <- ncol(design_matrix)
  if (any(cs != 0)) {
    error <- TRUE
    cs <- which(cs != 0)
    cs <- cat('Columns ', cs, 'do not sum to 0.\n')
  } else {
    cs <- ''
  }
  if (any(acs != 2)) {
    error <- TRUE
    acs <- which(acs != 2)
    acs <- cat('Columns ', acs, 'absolute value do not sum to 2.\n')
  } else {
    acs <- ''
  }
  if (cr != dn) {
    error <- TRUE
    cd <- cat(
      'Columns in design matrix (', dn,
      ') does not match the rows in contrast matrix (', cr, ')\n',
      'They need to match.')
  } else {
    cd <- ''
  }
  if (error) {
    stop(
      cat(cs, acs, cd)
    )
  }
}
