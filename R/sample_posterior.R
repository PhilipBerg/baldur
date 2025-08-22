utils::globalVariables(c("alpha", "betau", "id", "tmp", "intu", "condi"))
#' Sample the Posterior of the data and decision model and generate point
#' estimates
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#'
#'   Function to sample the posterior of the Bayesian data and
#'   decision model. It first produces the needed inputs for Stan's [sampling()]
#'   for each peptide (or protein, PTM, etc.). It then runs the sampling for the
#'   data and decision model. From the posterior, it then collects estimates and
#'   sampling statistics from the posterior of data model and integrates the
#'   decision distribution D. It then returns a [tibble()] with all the
#'   information for each peptide's posterior (see **Value** below). There are
#'   major time gains to be made by running this procedure in parallel.
#'   `infer_data_and_decision_model` has an efficient wrapper around
#'   `multidplyr`. This will let you to evenly distribute all peptides evenly to
#'   each worker. E.g., two workers will each run half of the peptides in
#'   parallel.
#'
#' @details The data model of Baldur assumes that the observations of a peptide,
#'   \eqn{\boldsymbol{Y}}, is a normally distributed with a standard deviation,
#'   \eqn{\sigma}, common to all measurements. In addition, it assumes that each
#'   measurement has a unique uncertainty \eqn{u}. It then models all
#'   measurements in the same condition with a common mean \eqn{\mu}. It then
#'   assumes that the common variation of the peptide is caused by the variation
#'   in the \eqn{\mu} As such, it models \eqn{\mu} with the common variance
#'   \eqn{\sigma} and a non-centered parametrization for condition level
#'   effects.
#' \deqn{
#'  \boldsymbol{Y}\sim\mathcal{N}(\boldsymbol{X}\boldsymbol{\mu},\sigma\boldsymbol{u})\quad
#'  \boldsymbol{\mu}\sim\mathcal{N}(\mu_0+\boldsymbol{\eta}\sigma,\sigma)
#' }
#'   It then assumes \eqn{\sigma} to be gamma distributed with hyperparameters
#'   infered from either a gamma regression [fit_gamma_regression] or a latent
#'   gamma mixture regression [fit_lgmr]. \deqn{\sigma\sim\Gamma(\alpha,\beta)}
#'   For details on the two priors for \eqn{\mu_0} see [empirical_bayes] or
#'   [weakly_informative]. Baldur then builds a posterior distribution of the
#'   difference(s) in means for contrasts of interest. In addition, Baldur
#'   increases the precision of the decision as the number of measurements
#'   increase. This is done by weighting the sample size with the contrast
#'   matrix. To this end, Baldur limits the possible contrast of interest such
#'   that each column must sum to zero, and the absolute value of each column
#'   must sum to two. That is, only mean comparisons are allowed.
#'   \deqn{
#'     \boldsymbol{D}\sim\mathcal{N}(\boldsymbol{\mu}^\text{T}\boldsymbol{K},\sigma\boldsymbol{\xi}),\quad \xi_{i}=\sqrt{\sum_{c=1}^{C}|k_{cm}|n_c^{-1}}
#'   }
#'   where \eqn{\boldsymbol{K}} is a contrast matrix of interest and
#'   \eqn{k_{cm}} is the contrast of the c:th condition in the m:th contrast of
#'   interest, and \eqn{n_c} is the number of measurements in the c:th
#'   condition. Baldur then integrates the tails of \eqn{\boldsymbol{D}} to
#'   determine the probability of error.
#'   \deqn{P(\text{\textbf{error}})=2\Phi(-\left|\boldsymbol{\mu}_{\boldsymbol{D}}-H_0\right|\odot\boldsymbol{\tau}_{\boldsymbol{D}})}
#'   where \eqn{H_0} is the null hypothesis for the difference in means,
#'   \eqn{\Phi} is the CDF of the standard normal,
#'   \eqn{\boldsymbol{\mu}_{\boldsymbol{D}}} is the posterior mean of
#'   \eqn{\boldsymbol{D}}, \eqn{\boldsymbol{\tau}_{\boldsymbol{D}}} is the
#'   posterior precision of \eqn{\boldsymbol{D}}, and \eqn{\odot} is the
#'   Hadamard product.
#'
#' @param data A `tibble` or `data.frame` with alpha,beta priors annotated
#' @param id_col A character for the name of the column containing the name of
#'   the features in data (e.g., peptides, proteins, etc.)
#' @param design_matrix A design matrix for the data. For the [empirical_bayes]
#'   prior only mean models are allowed (see example). For the
#'   [weakly_informative] prior more complicated design can be used.
#' @param contrast_matrix A contrast matrix of the decisions to test. Columns
#'   should sum to `0` and only mean comparisons are allowed. That is, the
#'   absolute value of the positive and negative values in each column has to
#'   sum to `2`. E.g., a column can be `[`0.5, 0.5, -1`]`\eqn{^T} but not `[`1,
#'   1, -1`]`\eqn{^T} or `[`0.5, 0.5, -2`]`\eqn{^T}. That is, `sum(abs(x))=2`
#'   where `x` is a column in the contrast matrix.
#' @param uncertainty_matrix A matrix of observation specific uncertainties
#' @param stan_model Which Bayesian model to use. Defaults to [empirical_bayes]
#'   but also allows [weakly_informative], or an user supplied function.
#' @param clusters The number of parallel threads/workers to run on.
#' @param h_not The value of the null hypothesis for the difference in means
#' @param auxiliary_columns Names of columns in the design matrix that does not
#' have corresponding data in the data set. For example, this can be
#' co-founding variables such as bashes.
#' @param ... Additional arguments passed to [sampling][rstan::sampling()].
#' Note that verbose will always be forced to `FALSE` to avoid console flooding.
#'
#' @return A [tibble][tibble::tibble()] or [data.frame()] annotated  with statistics of the
#'   posterior and p(error). `err` is the probability of error, i.e., the two
#'   tail-density supporting the null-hypothesis, `lfc` is the estimated
#'   \eqn{\log_2}-fold change, `sigma` is the common variance, and `lp` is the
#'   log-posterior. Columns without suffix shows the mean estimate from the
#'   posterior, while the suffixes `_025`, `_50`, and `_975`, are the 2.5, 50.0,
#'   and 97.5, percentiles, respectively. The suffixes `_eff` and `_rhat` are
#'   the diagnostic variables returned by `Stan`. In general, a larger `_eff`
#'   indicates effective sample size and `_rhat` indicates the mixing within
#'   chains and between the chains and should be smaller than 1.05 (please see
#'   the Stan manual for more details). In addition it returns a column
#'   `warnings` with potential warnings given by the sampler.
#'
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
#' contrast <- matrix(c(-1, 1), 2)
#' \donttest{
#' yeast_norm %>%
#'   head() %>% # Just running a few for the example
#'   infer_data_and_decision_model(
#'     'identifier',
#'     design,
#'     contrast,
#'     unc,
#'     clusters = 1 # I highly recommend increasing the number of parallel workers/clusters
#'                  # this will greatly reduce the time of running baldur
#'   )
#' }
infer_data_and_decision_model <- function(data,
                                          id_col,
                                          design_matrix,
                                          contrast_matrix,
                                          uncertainty_matrix,
                                          stan_model = empirical_bayes,
                                          clusters = 1,
                                          h_not = 0,
                                          auxiliary_columns = c(),
                                          ...) {

  rlang::is_missing(data)
  rlang::is_missing(id_col)
  rlang::is_missing(uncertainty_matrix)
  check_contrast_and_design(contrast_matrix, design_matrix, data, auxiliary_columns)
  check_id_col(id_col, colnames(data))

  rstan_inputs <- rlang::dots_list(...)
  rstan_inputs$verbose <- F

  N <- nrow(design_matrix)
  K <- ncol(design_matrix)
  C <- ncol(contrast_matrix)
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
                               "auxiliary_columns"
                             )
    )

    multidplyr::cluster_assign_partition(cl,
                                         id = ori_data[[id_col]],
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
                                                  auxiliary_columns
                                                )
                             )
    )
    results <- multidplyr::cluster_call(cl,
                                        purrr::map(dat,
                                                   purrr::quietly(stan_nse_wrapper),
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
                          rstan_inputs, auxiliary_columns
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

generate_stan_data_input <- function(id, id_col, design_matrix, data, uncertainty, comparison, N, K, C, alpha, beta, condi, condi_regex, auxiliary_columns){
  row <- data[data[[id_col]] == id, stringr::str_detect(names(data), condi_regex)]
  ybar <- purrr::map_dbl(purrr::set_names(condi, condi),
                         ~as.numeric(row[stringr::str_which(colnames(row), .x)]) %>%
                           mean()
  )[condi]
  ybar[auxiliary_columns] <- 0
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

stan_summary <- function(samp, dat, condi, contrast,  h_not){
  fit <- samp$result
  lfc_pars <- paste0("y_diff[", seq_len(ncol(contrast)), "]")
  sigma_pars <- paste0("sigma_lfc[", seq_len(ncol(contrast)), "]")
  err <- rstan::extract(fit, pars = lfc_pars) %>%
    purrr::map_dbl(estimate_error, h_not)
  summ <- rstan::summary(fit, pars = c("y_diff", 'sigma_lfc', 'lp__')) %>%
    magrittr::use_series(summary)
  summ <- summ[
    , colnames(summ) %in% c('mean', '2.5%', '50%', '97.5%', 'n_eff', 'Rhat')
  ]
  comps <- contrast_to_comps(contrast, condi)
  dplyr::tibble(
    comparison = comps,
    err = err,
    lfc = summ[rownames(summ) %in% lfc_pars, 1],
    lfc_025 = summ[rownames(summ) %in% lfc_pars, 2],
    lfc_50 = summ[rownames(summ) %in% lfc_pars, 3],
    lfc_975 = summ[rownames(summ) %in% lfc_pars, 4],
    lfc_eff = summ[rownames(summ) %in% lfc_pars, 5],
    lfc_rhat = summ[rownames(summ) %in% lfc_pars, 6],
    sigma = summ[rownames(summ) %in% sigma_pars, 1],
    sigma_025 = summ[rownames(summ) %in% sigma_pars, 2],
    sigma_50 = summ[rownames(summ) %in% sigma_pars, 3],
    sigma_975 = summ[rownames(summ) %in% sigma_pars, 4],
    sigma_eff = summ[rownames(summ) %in% sigma_pars, 5],
    sigma_rhat = summ[rownames(summ) %in% sigma_pars, 6],
    lp = summ[rownames(summ) == 'lp__', 1],
    lp_025 = summ[rownames(summ) == 'lp__', 2],
    lp_50 = summ[rownames(summ) == 'lp__', 3],
    lp_975 = summ[rownames(summ) == 'lp__', 4],
    lp_eff = summ[rownames(summ) == 'lp__', 5],
    lp_rhat = summ[rownames(summ) == 'lp__', 6],
    warnings = list(samp$warnings)
  )
}

bayesian_testing <- function(id, alpha, beta, data, id_col, design_matrix, comparison, model, N, K, C, uncertainty = NULL, h_not, rstan_args, auxiliary_columns){
  condi <- colnames(design_matrix)
  condi_regex <- get_conditions(design_matrix)
  dat <- generate_stan_data_input(id, id_col, design_matrix, data, uncertainty, comparison, N, K, C, alpha, beta, condi, condi_regex, auxiliary_columns)
  stan_output <- purrr::quietly(stan_nse_wrapper)(dat, model, rstan_args)
  stan_output %>%
    stan_summary(dat, condi, dat$c, h_not)
}

contrast_to_comps <- function(contrast, conditions){
  positives <- apply(contrast, 2, \(x) which(x>0)) %>%
    purrr::map_chr(~ stringr::str_flatten(conditions[.x], ' and '))
  negatives <- apply(contrast, 2, \(x) which(x<0)) %>%
    purrr::map_chr(~ stringr::str_flatten(conditions[.x], ' and '))
  stringr::str_c(positives, negatives, sep = ' vs ')
}

check_contrast_and_design <- function(contrast_matrix, design_matrix, data, auxiliary_columns) {
  cenv <- rlang::caller_call()
  design_matrix_no_aux <- check_design_aux(design_matrix, auxiliary_columns)
  check_contrast(contrast_matrix, cenv)
  check_design(design_matrix_no_aux, data, cenv)

  cr <- nrow(contrast_matrix)
  dn <- ncol(design_matrix)

  if (cr != dn) {
    rlang::abort(
      c(
        paste0(
          'There are ', dn, ' columns in design matrix but only ', cr,
          ' rows in contrast matrix'
        ),
        'They need to match.'),
      call = cenv
    )
  }
}

check_contrast <- function(contrast_matrix, cenv = rlang::caller_call()) {
  inp <- rlang::enquo(contrast_matrix)
  rlang::is_missing(contrast_matrix)
  check_matrix(inp, "contrast matrix", cenv)

  cs  <- !purrr::map_lgl(colSums(contrast_matrix),      ~ isTRUE(all.equal(.x, 0, check.attributes = F)))
  acs <- !purrr::map_lgl(colSums(abs(contrast_matrix)), ~ isTRUE(all.equal(.x, 2, check.attributes = F)))

  if (any(cs)) {
    csi <- which(cs)
    rlang::abort(
      c(
        paste0(
          "Your contrast matix, ",
          rlang::as_name(inp),
          ", column(s) do not sum to zero.")
        ,
        paste0(
          'The following column(s): ',
          stringr::str_flatten(csi, ', '),
          ', do not sum to 0.'
        ),
        paste0(
          'They sum to ', stringr::str_flatten(cs[csi], ', '), ', respesctively.'
        )
      ),
      call = cenv
    )
  }
  if (any(acs)) {
    acsw <- which(acs)
    rlang::abort(
      c(
        paste0('Contrast column(s) ', acsw, "s absolute value do not sum to 2."),
        paste0('They sum to ', acs[acsw])
      ),
      call = cenv
    )
  }
}
