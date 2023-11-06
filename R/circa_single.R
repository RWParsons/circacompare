#' @title circa_single
#' @name circa_single
#'
#' @description \code{circa_single} performs an analysis on a single rhythmic dataset. It estimates the mesor, amplitude and phase of the data provided.
#'
#' @param x data.frame.  This is the data.frame which contains the rhythmic data in a tidy format.
#' @param col_time The name of the column within the data.frame, x, which contains time in hours at which the data were collected.
#' @param col_outcome The name of the column within the data.frame, x, which contains outcome measure of interest.
#' @param period The period of the rhythm. For circadian rhythms, leave this as the default value, 24.
#' @param alpha_threshold The level of alpha for which the presence of rhythmicity is considered. Default is 0.05.
#' @param timeout_n The upper limit for the model fitting attempts. Default is 10,000.
#' @param return_figure Whether or not to return a ggplot graph of the rhythm and cosine model.
#' @param control \code{list}. Used to control the parameterization of the model.
#' @param weights An optional numeric vector of (fixed) weights. When present, the objective function is weighted least squares.
#' @param suppress_all Logical. Set to \code{TRUE} to avoid seeing errors or messages during model fitting procedure. Default is \code{FALSE}.
#'
#' @return list
#' @export
#'
#' @examples
#' df <- make_data()
#' df <- df[df$group == "g1", ]
#' out <- circa_single(x = df, col_time = "time", col_outcome = "measure")
#' out
#'
#' # with sample weights (arbitrary weights for demonstration)
#' sw <- runif(n = nrow(df))
#' out2 <- circa_single(
#'   x = df,
#'   col_time = "time",
#'   col_outcome = "measure",
#'   weights = sw,
#'   suppress_all = TRUE
#' )
#' out2
#'
circa_single <- function(x,
                         col_time,
                         col_outcome,
                         period = 24,
                         alpha_threshold = 0.05,
                         timeout_n = 10000,
                         return_figure = TRUE,
                         control = list(),
                         weights = NULL,
                         suppress_all = FALSE) {
  controlVals <- circa_single_control()
  controlVals[names(control)] <- control

  if (controlVals$period_param) {
    controlVals$main_params <- c(controlVals$main_params, "tau")
  }
  if ("tau" %in% controlVals$main_params) {
    controlVals$period_param <- TRUE
  }

  x <- x[c(col_time, col_outcome)]
  colnames(x) <- c("time", "measure")

  if (!class(x$time) %in% c("numeric", "integer")) {
    stop(paste("The time variable which you gave was a '",
      class(x$time), "' \nThis function expects time to be given as hours and be of class 'integer' or 'numeric'.",
      "\nPlease convert the time variable in your dataframe to be of one of these classes",
      sep = ""
    ))
  }

  if (!class(x$measure) %in% c("numeric", "integer")) {
    stop(paste("The measure variable which you gave was a '",
      class(x$measure), "' \nThis function expects measure to be number and be of class 'integer' or 'numeric'.",
      "\nPlease convert the measure variable in your dataframe to be of one of these classes",
      sep = ""
    ))
  }

  if (!class(period) %in% c("numeric", "integer") & !controlVals$period_param) {
    stop(paste0(
      "The period argument must be a number representing the period of the rhythm in hours\n",
      "If you would like the period to be estimated as part of the model, use:\ncontrol=list(period_param=TRUE)"
    ))
  }

  if (controlVals$period_param & !is.na(period)) {
    message(paste0("control$period_param is TRUE\n'period=", period, "' is being ignored.\nSet 'period=NA' to avoid this message"))
  }
  x$time_r <- x$time * 2 * pi
  if (!controlVals$period_param) {
    x$period <- period
  } else {
    if (is.null(controlVals$period_min) | is.null(controlVals$period_min)) {
      message(paste0(
        "If you want the model to estimate the period using a parameter,",
        "you may get faster convergence if you provide an approximate range using 'period_min' and 'period_max' in control()",
        "\nCurrently assuming period is between: period_min=", controlVals$period_min,
        "and period_max=", controlVals$period_max
      ))
    }
  }

  if (!is.null(weights)) {
    check_weights(x, weights)
    x$weights <- weights
  } else {
    x$weights <- rep(1, nrow(x))
  }

  success <- FALSE
  n <- 0
  form <- create_formula(main_params = controlVals$main_params, decay_params = controlVals$decay_params)$formula

  while (!success) {
    fit.nls <- try(
      {
        stats::nls(
          formula = form,
          data = x,
          start = start_list(outcome = x$measure, controlVals = controlVals),
          weights = weights
        )
      },
      silent = suppress_all
    )
    if (inherits(fit.nls, "try-error")) {
      n <- n + 1
    } else {
      nls_coefs <- extract_model_coefs(fit.nls)
      V <- nls_coefs[, "estimate"]
      success <- assess_model_estimates(param_estimates = V, controlVals = controlVals)
      n <- n + 1
    }
    if (n >= timeout_n) {
      stop("Failed to converge data prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript).")
    }
  }
  data_rhythmic <- nls_coefs["alpha", "p_value"] < alpha_threshold

  if (!controlVals$period_param) {
    V["tau"] <- period
  }

  eq_expression <- create_formula(main_params = controlVals$main_params, decay_params = controlVals$decay_params)$f_equation
  eval(parse(text = eq_expression))


  if (return_figure) {
    p <- ggplot2::ggplot(x, ggplot2::aes(time, measure)) +
      ggplot2::geom_point() +
      ggplot2::xlim(
        min(floor(x$time / period) * period),
        max(ceiling(x$time / period) * period)
      )
    if (data_rhythmic) {
      fig_out <- p +
        ggplot2::stat_function(fun = eq, size = 1) +
        ggplot2::labs(subtitle = "Data is rhythmic", x = "time (hours)")
    } else {
      fig_out <- p +
        ggplot2::labs(subtitle = "Data is arrhythmic", x = "time (hours)")
    }
  }

  results_summary <-
    circa_summary(model = fit.nls, period = period, control = controlVals)

  if (return_figure) {
    return(list(fit = fit.nls, summary = results_summary, plot = fig_out))
  } else {
    return(list(fit = fit.nls, summary = results_summary))
  }
}


circa_single_control <- function(period_param = F, period_min = 20, period_max = 28,
                                 main_params = c("k", "alpha", "phi"), decay_params = c()) {
  list(
    period_param = period_param, period_min = period_min, period_max = period_max,
    main_params = main_params, decay_params = decay_params
  )
}
