#' Fit Latent Gamma Mixture Regression
#'
#' @param data A `data.frame` with mean-variance trends to use in the fitting.
#'  The columns need to have the following hard-coded names: `mean` and `sd`.
#' @param model Defaults to [lgmr_model], can also be user supplied
#'   stan_model
#' @param iter Total number of samples to draw
#' @param warmup Number of warm-up samples to draw
#' @param chains Number of chains to run
#' @param cores Number of cores to use per chain
#' @param return_stanmodel Should the `stanfit` object be returned with the
#'   model?
#' @param simplify Should only the mean estimates of the posterior be returned?
#' @param ... Additional arguments to `rstan`'s [sampling][rstan::sampling()].
#'
#' @return A fitted `lgmr` model.
#' @export
#' @name fit_lgmr
#'
#' @examples
#' # Define design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#'\donttest{
#' # Normalize data, calculate M-V trend, and fit LGMR model
#' yeast_lgmr <- yeast %>%
#'     psrn("identifier") %>%
#'     # Get mean-variance trends
#'     calculate_mean_sd_trends(design) %>%
#'     # Fit the model
#'     fit_lgmr(lgmr_model)
#' }
fit_lgmr <- function(data, model, iter = 6000, warmup = 1500, chains = 5, cores = 1, return_stanmodel = FALSE, simplify = FALSE, ...) {
  if (!"sd" %in% names(data)) {
    stop("sd is not a column in the data\n Did you forget to calculate the Mean-Variance trend?")
  } else if(is.unsorted(data$sd) & is.unsorted(data$sd, strictly = TRUE)) {
    stop("Data is not sorted from smallest to largest sd.\n Did you forget to sort the data?")
  }
  stan_args <- list(
    control = list(
      adapt_delta = .9
    )
  )

  input_args <- rlang::dots_list()
  stan_args[names(input_args)] <- input_args
  n <- nrow(data)
  input <- list(
    N = n,
    y = data$sd, x = data$mean
  )

  samp <- rstan::sampling(
    model, data = input,
    iter = iter, warmup = warmup,
    chains = chains, cores = cores
  )

  summ <- rstan::summary(samp) %>%
    use_series(summary)

  fitted_model       <- structure(list(), class = "lgmr")
  fitted_model$coef  <- model_summary[rownames(model_summary) %in% c("I", "S", "I_L", "S_L"),]
  fitted_model$aux   <- model_summary[rownames(model_summary) %in% c("alpha", "nrmse"),]
  fitted_model$theta <- model_summary[stringr::str_detect(rownames(model_summary), "theta"),]
  if (simplify) {
    fitted_model <- purrr::map(fitted_model, ~ .x[,"mean"])
  }
  if (return_stanmodel) {
    fitted_model$stanfit <- samp
  }
  return(fitted_model)
}


#' @param x A `lgrm` model.
#' @param pars If you want to print the regression coefficients, theta, or both.
#' @param digits Number digits to print
#'
#' @rdname fit_lgmr
#' @export
print.lgmr <- function(x, pars = c("coefficients", "theta"), digits = 3) {
  mu <- round(x$coef[, 'mean'], digits)
  cat("\nLGMR Model\n")
  pars <- match.arg(pars, several.ok = TRUE)
  cat("\tmu=", "exp(", mu["I"],
      " - ",
      mu["S"], " f(bar_y)) + kappa exp(",
      mu["I_L"],
      " - ",
      mu["S_L"],
      " f(bar_y))", sep = ''
  )
  if (pars == "coefficients") {
    cat("\n\n",
        "Coefficients:\n"
    )
    print.default(
      x$reg,
      digits = digits,
      print.gap = 2L,
      quote = FALSE
    )
    cat("\n\n",
        "Auxilary:\n"
    )
    print.default(
      x$aux,
      digits = digits,
      print.gap = 2L,
      quote = FALSE
    )
  }
  if (pars == "theta") {
    cat("\n\n",
        "theta:\n"
    )
    print.default(
      x$theta,
      digits = digits,
      print.gap = 2L,
      quote = FALSE
    )
  }
}
