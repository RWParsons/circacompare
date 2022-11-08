#' @title circacompare_mixed
#' @name circacompare_mixed
#'
#' @description \code{circacompare_mixed} is similar to \code{circacompare} but allows for some simple, user-specified random-effects on the rhythmic parameters of choice.
#'
#' @param x \code{data.frame}.  This is the data.frame which contains the rhythmic data for two groups in a tidy format.
#' @param col_time The name of the column within the data.frame, \code{x}, which contains time in hours at which the data were collected.
#' @param col_group The name of the column within the data.frame, \code{x}, which contains the grouping variable.  This should only have two levels.
#' @param col_outcome The name of the column within the data.frame, \code{x}, which contains outcome measure of interest.
#' @param col_id The name of the column within the data.frame, \code{x}, which contains the identifying values for the random effect, such as \code{subject_id}.
#' @param randomeffects which rhythmic parameters to allow random effects. The default is to include no rhythmic parameters.
#' @param period The period of the rhythm. For circadian rhythms, leave this as the default value, \code{24}.
#' @param alpha_threshold The level of alpha for which the presence of rhythmicity is considered. Default is to \code{0.05}.
#' @param nlme_control A list of control values for the estimation algorithm to replace the default values returned by the function nlme::nlmeControl. Defaults to an empty list.
#' @param nlme_method A character string. If "REML" the model is fit by maximizing the restricted log-likelihood. If "ML" the log-likelihood is maximized. Defaults to "REML".
#' @param suppress_all Logical. Set to \code{TRUE} to avoid seeing errors or messages during model fitting procedure. Default is \code{FALSE}. If \code{FALSE}, also runs \code{nlme()} with \code{verbose = TRUE}.
#' @param timeout_n The upper limit for the model fitting attempts. Default is \code{10000}.
#' @param control \code{list}. Used to control the parameterization of the model.
#'
#' @return list
#' @export
#'
#' @examples
#' # Generate some data with within-id correlation for phase-shift (phi1)
#' set.seed(99)
#' phi1_in <- 3.15
#'
#' mixed_data <- function(n) {
#'   counter <- 1
#'   for (i in 1:n) {
#'     x <- make_data(k1 = 0, alpha1 = 0, phi1 = rnorm(1, phi1_in, 0.5), hours = 72, noise_sd = 1)
#'     x$id <- counter
#'     counter <- counter + 1
#'     if (i == 1) {
#'       res <- x
#'     } else {
#'       res <- rbind(res, x)
#'     }
#'   }
#'   return(res)
#' }
#' df <- mixed_data(20)
#' out <- circacompare_mixed(
#'   x = df,
#'   col_time = "time",
#'   col_group = "group",
#'   col_outcome = "measure",
#'   col_id = "id",
#'   control = list(grouped_params = c("phi"), random_params = c("phi1"))
#' )
#' out
#'
circacompare_mixed <- function(x,
                               col_time,
                               col_group,
                               col_outcome,
                               col_id,
                               randomeffects = c(),
                               period = 24,
                               alpha_threshold = 0.05,
                               nlme_control = list(),
                               nlme_method = "REML",
                               suppress_all = FALSE,
                               timeout_n = 10000,
                               control = list()) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("{nlme} is required for circacompare_mixed(), please install {nlme}")
  }

  controlVals <- circacompare_mixed_control()
  controlVals[names(control)] <- control

  if (controlVals$period_param & !"tau" %in% controlVals$main_params) {
    controlVals$main_params <- c(controlVals$main_params, "tau")
  }
  if ("tau" %in% controlVals$main_params) {
    controlVals$period_param <- TRUE
  }

  randomeffects <- controlVals$random_params <- unique(c(randomeffects, controlVals$random_params))

  x <- x[c(col_time, col_group, col_id, col_outcome)]
  colnames(x) <- c("time", "group", "id", "measure")

  if (length(controlVals$decay_params) > 0) {
    p <- append(controlVals$main_params, paste0(controlVals$decay_params, "_decay"))
  } else {
    p <- controlVals$main_params
  }
  controlVals$non_grouped_params <- p
  p <- append(p, paste0(controlVals$grouped_params, "1"))

  if (length(setdiff(randomeffects, p)) != 0) {
    stop('"randomeffects" should only include the names of parameters\nthat represent rhythmic characteristics in the model.\nThey should be a subset of, or equal to c("k", "k1", "alpha", "alpha1", "phi", "phi1")')
  }

  if (length(randomeffects) == 0) {
    stop("If you do not want to include any random effects, than you ought to use 'circacompare' rather than 'circacompare_mixed'")
  }

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
  randomeffects_single <- intersect(controlVals$non_grouped_params, controlVals$random_params)

  g1_model <- model_each_group(
    data = dat_group_1, type = controlVals$fit_group_models_method, form = form_single,
    controlVals = controlVals,
    args = list(
      timeout_n = timeout_n,
      alpha_threshold = alpha_threshold,
      randomeffects = randomeffects_single,
      nlme_method = nlme_method,
      nlme_control = nlme_control,
      suppress_all = suppress_all
    )
  )

  if (g1_model$timeout) {
    stop("Failed to converge", group_1_text, " model prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript).")
  }

  g2_model <- model_each_group(
    data = dat_group_2, type = controlVals$fit_group_models_method, form = form_single,
    controlVals = controlVals,
    args = list(
      timeout_n = timeout_n,
      alpha_threshold = alpha_threshold,
      randomeffects = randomeffects_single,
      nlme_method = nlme_method,
      nlme_control = nlme_control,
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

  n <- 0
  success <- FALSE
  form_group <- create_formula(main_params = controlVals$main_params, decay_params = controlVals$decay_params, grouped_params = controlVals$grouped_params)$formula
  fixedeffects_formula <- stats::formula(paste(paste0(p, collapse = "+"), "~ 1"))
  randomeffects_formula <- stats::formula(paste(paste0(randomeffects, collapse = "+"), "~ 1 | id"))

  while (!success) {
    fit.nlme <- try(
      {
        nlme::nlme(
          model = form_group,
          random = randomeffects_formula,
          fixed = fixedeffects_formula,
          data = x,
          start = unlist(start_list_grouped(g1 = g1_model$model, g2 = g2_model$model, grouped_params = controlVals$grouped_params)),
          method = nlme_method,
          control = nlme_control,
          verbose = !suppress_all
        )
      },
      silent = ifelse(suppress_all, TRUE, FALSE)
    )

    if ("try-error" %in% class(fit.nlme)) {
      n <- n + 1
    } else {
      nlme_coefs <- extract_model_coefs(fit.nlme)
      V <- nlme_coefs[, "estimate"]
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
      min(floor(x$time / period) * period),
      max(ceiling(x$time / period) * period)
    )

  results_summary <-
    circa_summary(
      model = fit.nlme, period = period, control = controlVals,
      g1 = g1_model, g2 = g2_model, g1_text = group_1_text, g2_text = group_2_text
    )

  return(list(plot = fig_out, summary = results_summary, fit = fit.nlme))
}

circacompare_mixed_control <- function(period_param = F, period_min = 20, period_max = 28,
                                       main_params = c("k", "alpha", "phi"),
                                       grouped_params = c("k", "alpha", "phi"),
                                       decay_params = c(), random_params = c(),
                                       fit_group_models_method = "nlme") {
  list(
    period_param = period_param, period_min = period_min, period_max = period_max,
    main_params = main_params, grouped_params = grouped_params,
    decay_params = decay_params, random_params = random_params,
    fit_group_models_method = fit_group_models_method
  )
}
