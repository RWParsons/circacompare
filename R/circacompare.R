#' @title circacompare
#' @name circacompare
#'
#' @description \code{circacompare} performs a comparison between two rhythmic groups of data. It tests for rhythmicity and then fits a nonlinear model with parametrization to estimate and statistically support differences in mesor, amplitude, and phase between groups.
#'
#' @param x data.frame.  This is the data.frame which contains the rhythmic data for two groups in a tidy format.
#' @param col_time The name of the column within the data.frame, x, which contains time in hours at which the data were collected.
#' @param col_group The name of the column within the data.frame, x, which contains the grouping variable.  This should only have two levels.
#' @param col_outcome The name of the column within the data.frame, x, which contains outcome measure of interest.
#' @param period The period of the rhythm. For circadian rhythms, leave this as the default value, 24.
#' @param alpha_threshold The level of alpha for which the presence of rhythmicity is considered. Default is 0.05.
#' @param timeout_n The upper limit for the model fitting attempts. Default is 10,000.
#' @param control \code{list}. Used to control the parameterization of the model.
#' @param sample_weights a numeric vector of per-sample weights (one per sample) to allow downweighting of outliers in the differential analysis
#' @param suppress_all Logical. Set to \code{TRUE} to avoid seeing errors or messages during model fitting procedure. Default is \code{FALSE}.
#'
#' @return list
#' @export
#'
#' @examples
#' df <- make_data(phi1 = 6)
#' out <- circacompare(
#'   x = df, col_time = "time", col_group = "group",
#'   col_outcome = "measure"
#' )
#' out
#' # with sample weights
#' set.seed(1)
#' sw <- jitter(rep(1, nrow(df)), factor=2)
#' out <- circacompare(
#'   x = df, col_time = "time", col_group = "group",
#'   col_outcome = "measure", sample_weights = sw
#' )
circacompare <- function(x,
                         col_time,
                         col_group,
                         col_outcome,
                         period = 24,
                         alpha_threshold = 0.05,
                         timeout_n = 10000,
                         control = list(),
                         sample_weights = NULL,
                         suppress_all = FALSE) {
  controlVals <- circacompare_control()
  controlVals[names(control)] <- control

  if (controlVals$period_param & !"tau" %in% controlVals$main_params) {
    controlVals$main_params <- c(controlVals$main_params, "tau")
  }
  if ("tau" %in% controlVals$main_params) {
    controlVals$period_param <- TRUE
  }

  x <- x[c(col_time, col_group, col_outcome)]
  colnames(x) <- c("time", "group", "measure")

  if (length(levels(as.factor(x$group))) != 2) {
    stop("Your grouping variable had more or less than 2 levels! \nThis function is used to compare two groups of data. \nTo avoid me having to guess, please send data with only two possible values in your grouping variable to this function.")
  }

  if (!class(x$time) %in% c("numeric", "integer")) {
    stop(paste("The time variable which you gave was a '",
               class(x$time),
               "' \nThis function expects time to be given as hours and be of class 'integer' or 'numeric'.",
               "\nPlease convert the time variable in your dataframe to be of one of these classes",
               sep = ""
    ))
  }

  if (!class(x$measure) %in% c("numeric", "integer")) {
    stop(paste("The measure variable which you gave was a '",
               class(x$measure),
               "' \nThis function expects measure to be number and be of class 'integer' or 'numeric'.",
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

  group_1_text <- levels(as.factor(x$group))[1]
  group_2_text <- levels(as.factor(x$group))[2]
  x$x_group <- ifelse(x$group == group_1_text, 0, 1)
  dat_group_1 <- x[x$group == group_1_text, ]
  dat_group_2 <- x[x$group == group_2_text, ]

  form_single <- create_formula(main_params = controlVals$main_params, decay_params = controlVals$decay_params)$formula

  g1_model <- model_each_group(
    data = dat_group_1, type = "nls", form = form_single,
    controlVals = controlVals,
    args = list(
      timeout_n = timeout_n,
      alpha_threshold = alpha_threshold,
      suppress_all = suppress_all
    )
  )

  if (g1_model$timeout) {
    stop("Failed to converge", group_1_text, " model prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript).")
  }

  g2_model <- model_each_group(
    data = dat_group_2, type = "nls", form = form_single,
    controlVals = controlVals,
    args = list(
      timeout_n = timeout_n,
      alpha_threshold = alpha_threshold,
      suppress_all = suppress_all
    )
  )
  if (g2_model$timeout) {
    stop("Failed to converge", group_2_text, " model prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript).")
  }

  both_groups_rhythmic <- ifelse(g1_model$rhythmic & g2_model$rhythmic, TRUE, FALSE)

  if (!both_groups_rhythmic) {
    if (!g1_model$rhythmic & !g2_model$rhythmic) {
      stop("Both groups of data were arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups.")
    }
    if (!g1_model$rhythmic) {
      stop(group_1_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups.")
    } else {
      stop(group_2_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups.")
    }
  }

  if(!is.null(sample_weights)){
    l.weights <- length(sample_weights)
    if(l.weights != nrow(x) | sum(is.na(sample_weights)) > 0 | sum(sample_weights <= 0) > 0 | !is.numeric(sample_weights))
      stop("sample_weights must be numeric, positive, non-zero, non-NA and of same length as the number of rows in x")
    x$sample_weights <- sample_weights
  } else x$sample_weights <- rep(1, nrow(x))

  n <- 0
  success <- FALSE
  form_group <- create_formula(main_params = controlVals$main_params, decay_params = controlVals$decay_params, grouped_params = controlVals$grouped_params)$formula
  while (!success) {
    fit.nls <- try(
      {
        stats::nls(
          formula = form_group,
          data = x,
          start = start_list_grouped(g1 = g1_model$model, g2 = g2_model$model, grouped_params = controlVals$grouped_params),
          control = stats::nls.control(maxiter = 100, minFactor = 1 / 10000),
          weights = sample_weights
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
    if (n > timeout_n) {
      stop("Both groups of data were rhythmic but the curve fitting procedure failed due to timing out. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript).")
    }
  }
  if (!controlVals$period_param) {
    V["tau"] <- period
  }

  eq_expression <- create_formula(
    main_params = controlVals$main_params,
    decay_params = controlVals$decay_params,
    grouped_params = controlVals$grouped_params
  )$f_equation

  eval(parse(text = eq_expression$g1))
  eval(parse(text = eq_expression$g2))

  fig_out <- ggplot2::ggplot(x, ggplot2::aes(time, measure)) +
    ggplot2::stat_function(ggplot2::aes(colour = group_1_text), fun = eq_1, size = 1) +
    ggplot2::stat_function(ggplot2::aes(colour = group_2_text), fun = eq_2, size = 1) +
    ggplot2::geom_point(ggplot2::aes(colour = group)) +
    ggplot2::scale_colour_manual(
      breaks = c(group_1_text, group_2_text),
      values = c("blue", "red")
    ) +
    ggplot2::labs(
      colour = "Legend",
      x = "time (hours)"
    ) +
    ggplot2::xlim(
      min(floor(x$time / V["tau"]) * V["tau"]),
      max(ceiling(x$time / V["tau"]) * V["tau"])
    )

  results_summary <-
    circa_summary(
      model = fit.nls, period = period, control = controlVals,
      g1 = g1_model, g2 = g2_model, g1_text = group_1_text, g2_text = group_2_text
    )

  return(list(plot = fig_out, summary = results_summary, fit = fit.nls))
}


circacompare_control <- function(period_param = F, period_min = 20, period_max = 28,
                                 main_params = c("k", "alpha", "phi"),
                                 grouped_params = c("k", "alpha", "phi"),
                                 decay_params = c()) {
  list(
    period_param = period_param, period_min = period_min, period_max = period_max,
    main_params = main_params, grouped_params = grouped_params,
    decay_params = decay_params
  )
}
