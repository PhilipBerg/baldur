utils::globalVariables(c("reg", "coef"))
#'Fit Latent Gamma Mixture Regression
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' See [lgmr_model] for model details.
#'
#'@param data A `data.frame` with mean-variance trends to use in the fitting.
#'  The columns need to have the following hard-coded names: `mean` and `sd`.
#'@param  id_col A character for the name of the column containing the name of
#'  the features in data (e.g., peptides, proteins, etc.)
#'@param model Defaults to [lgmr_model] (see it for details on the model), can
#'  also be an user supplied [stan_model][rstan::stan_model()]
#'@param iter Total number of samples to draw
#'@param warmup Number of warm-up samples to draw
#'@param chains Number of chains to run
#'@param cores Number of cores to use per chain
#'@param return_stanfit Should the `stanfit` object be returned with the model?
#'@param simplify Should only the mean estimates of the posterior be returned?
#'@param ... Additional arguments to `rstan`'s [sampling][rstan::sampling()].
#'  Does nothing for `print` or `coef` only for `fit_lgmr`.
#'@param id_col A character for the name of the column containing the name of
#'  the features in data (e.g., peptides, proteins, etc.). Has to be a unique
#'  identifier for each feature.
#'
#'@return A fitted `lgmr` model.
#'@export
#'@name fit_lgmr
#'
#' @examples
#' # Define design matrix
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
#'
#' }
fit_lgmr <- function(data, id_col, model = lgmr_model, iter = 6000, warmup = 1500, chains = 5, cores = 1, return_stanfit = FALSE, simplify = FALSE, ...) {
  if (!"sd" %in% names(data)) {
    stop("sd is not a column in the data\n Did you forget to calculate the Mean-Variance trend?")
  } else if(!"mean" %in% names(data)) {
    stop("mean is not a column in the data\n Did you forget to calculate the Mean-Variance trend?")
  } else if (missing(id_col)) {
    stop('id_col is a required input.')
  }
  check_id_col(id_col, colnames(data))

  stan_args <- list(
    control = list(
      adapt_delta = .9
    )
  )

  n <- nrow(data)
  input <- list(
    N = n,
    y = data$sd, x = data$mean
  )

  input_args <- rlang::dots_list(...)
  stan_args[names(input_args)] <- input_args
  samp <- rlang::eval_tidy(
    rlang::call2(rstan::sampling,
                 model, data = input,
                 iter = iter, warmup = warmup,
                 chains = chains, cores = cores,
                 !!!stan_args
    )
  )

  model_summary <- rstan::summary(samp) %>%
    use_series(summary)

  fitted_model                 <- list()
  fitted_model$coef            <- model_summary[rownames(model_summary) %in% c("I", "S", "I_L", "S_L"),]
  fitted_model$aux             <- model_summary[rownames(model_summary) %in% c("alpha", "nrmse"),]
  fitted_model$theta           <- model_summary[stringr::str_detect(rownames(model_summary), "theta"),]
  rownames(fitted_model$theta) <- paste0('theta_', data[[id_col]])

  if (simplify) {
    fitted_model <- purrr::map(fitted_model, ~ .x[,"mean"])
  }

  fitted_model$stan_model <- model
  fitted_model$data       <- dplyr::select(data, mean, sd)
  fitted_model$simplify   <- simplify

  if (return_stanfit) {
    fitted_model$stanfit <- samp
  }

  fitted_model            <- structure(fitted_model, class = "lgmr")

  return(fitted_model)
}


#' @param x,object An `lgmr` model.
#' @param pars If you want to print/extract the regression coefficients, theta, auxiliary (alpha and NRMSE), or all
#' @param digits Number of digits to print
#'
#' @rdname fit_lgmr
#' @export
print.lgmr <- function(x, simplify = x$simplify, pars = c("auxiliary", "coefficients"), digits = 3, ...) {

  if (!simplify & x$simplify) {
    options(warn = 1)
    rlang::warn(
      c(
        "Model is already simplified, cannot print full model specs.",
        "Defaulting to simplify = TRUE"
      )
    )
    options(warn = 0)
    simplify <- TRUE
  }

  mu <- round(coef(x, TRUE, "coefficients"), digits = digits)
  pars <- match_pars(pars)
  cf <- coef(x, simplify, pars)
  cat("\nLGMR Model\n")
  cat("\t\U03bc = ", "exp(", mu["I"],
      " - ",
      mu["S"], " f(\U0079\U0304)) + \U03ba exp(\U03b8(",
      mu["I_L"],
      " - ",
      mu["S_L"],
      " f(\U0079\U0304)))", sep = ''
  )
  if (is.list(cf)) {
    cf %>%
      purrr::imap(set_print_names, simplify) %>%
      purrr::walk2(pars, build_print, digits)
  } else {
    set_print_names(cf, pars, simplify) %>%
      build_print(pars, digits)
  }
}

set_print_names <- function(x, type, simplify) {
  replacement <- dplyr::case_when(
    type == "coef" ~ list(setNames(
      c("\U03b3_0L", "\U03b3_0", "\U03b3_\U0079\U0304L", "\U03b3_\U0079\U0304"),
      c("I_L", "I", "S_L", "S")
    )),
    type == "aux" ~ list(setNames(c("\U03b1", "NRMSE"), c('alpha', 'nrmse'))),
    T             ~ list(setNames("\U03b8", "theta"))
  )[[1]]
  n <- list()
  if (simplify) {
    n$get <-  names
    n$set <- `names<-`
  } else {
    n$get <-  rownames
    n$set <- `rownames<-`
  }
  x <- n$set(x,
        stringr::str_replace_all(
          n$get(x),
          replacement
        )
  )
  return(x)
}

build_print <- function(x, type, digits) {
  top <- dplyr::case_match(
    type,
    "auxiliary"    ~ "Auxiliary:\n",
    "coefficients" ~ "Coefficients:\n",
    .default       = "Latent contribution (\U03b8):\n"
  )
  cat("\n\n", top)
  print.default(
    x,
    digits = digits,
    print.gap = 2L,
    quote = FALSE
  )
}

#' @rdname fit_lgmr
#'
#' @export
coef.lgmr <- function(object, simplify = FALSE, pars = c("coefficients", "auxiliary"), ...) {

  pars <- match_pars(pars)
  if (simplify & !object$simplify) {
    f <- function(x) x[, "mean"]
  } else {
    f <- function(x) x
  }

  vars <- dplyr::case_when(
    pars == "coefficients" ~ 'coef',
    pars == "auxiliary" ~ 'aux',
    pars == "theta" ~ 'theta'
  ) %>%
    setNames(., .) %>%
    map(
      ~ f(object[[.x]])
    )

  if (length(vars) == 1) {
    vars <- vars[[1]]
  }

  return(vars)
}

mu_fun <- function(theta, reg_pars, y_bar, m, s){
  check_mu_fun_inputs(theta, y_bar)
  y_bar_star <- (y_bar - m)/s
  exp(reg_pars['I'] - reg_pars['S']*y_bar_star) +
    0.001 * exp(theta * (reg_pars["I_L"] - reg_pars["S_L"] * y_bar_star))
}

match_pars <- function(pars = c("coefficients", "auxiliary", "theta", "all")) {
  if ("all" %in% pars) {
    c("auxiliary", "coefficients", "theta")
  } else {
    match.arg(pars, several.ok = TRUE)
  }
}

check_mu_fun_inputs <- function(theta, y_bar) {
  if (length(theta) != length(y_bar)) {
    rlang::abort(
      c(
        "Length of theta not equal to length of y_bar.",
        "mu_fun does not recycle.",
        "If you did not run mu_fun yourself; please report this error to: https://github.com/PhilipBerg/baldur/issues"
      )
    )
  }
}
