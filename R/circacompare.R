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
#'
#' @return list
#' @export
#'
#' @examples
#' df <- make_data(phi1 = 6)
#' out <- circacompare(x = df, col_time = "time", col_group = "group",
#'                     col_outcome = "measure")
#' out
circacompare <- function(x,
                         col_time,
                         col_group,
                         col_outcome,
                         period = 24,
                         alpha_threshold = 0.05,
                         timeout_n = 10000){
  if(!"ggplot2" %in% utils::installed.packages()[, "Package"]){
    return(message("Please install 'ggplot2'"))
  }

  x <- x[c(col_time, col_group, col_outcome)]
  colnames(x) <- c("time", "group", "measure")

  if(length(levels(as.factor(x$group))) != 2){
    return(message("Your grouping variable had more or less than 2 levels! \nThis function is used to compare two groups of data. \nTo avoid me having to guess, please send data with only two possible values in your grouping variable to this function."))
  }

  group_1_text <- levels(as.factor(x$group))[1]
  group_2_text <- levels(as.factor(x$group))[2]

  if(!class(x$time) %in% c("numeric", "integer")){
    return(message(paste("The time variable which you gave was a '",
                         class(x$time),
                         "' \nThis function expects time to be given as hours and be of class 'integer' or 'numeric'.",
                         "\nPlease convert the time variable in your dataframe to be of one of these classes",
                         sep = "")))
  }

  if(!class(x$measure) %in% c("numeric", "integer")){
    return(message(paste("The measure variable which you gave was a '",
                         class(x$measure),
                         "' \nThis function expects measure to be number and be of class 'integer' or 'numeric'.",
                         "\nPlease convert the measure variable in your dataframe to be of one of these classes",
                         sep = "")))
  }

  x$time_r <- (x$time/24)*2*pi*(24/period)
  x$x_group <- ifelse(x$group == group_1_text, 0, 1)

  comparison_model_success <- 0
  comparison_model_timeout <- FALSE
  g1_success <- 0
  g2_success <- 0
  g1_alpha_p <- NA
  g2_alpha_p <- NA
  n <- 0
  dat_group_1 <- x[x$group == group_1_text,]
  dat_group_2 <- x[x$group == group_2_text,]

  g1_model <- model_each_group(data=dat_group_1, type="nls",
                               args=list(
                                 timeout_n=timeout_n,
                                 alpha_threshold=alpha_threshold
                               ))
  if(g1_model$timeout){return(message("Failed to converge", group_1_text, " model prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))}

  g2_model <- model_each_group(data=dat_group_1, type="nls",
                               args=list(
                                 timeout_n=timeout_n,
                                 alpha_threshold=alpha_threshold
                               ))
  if(g2_model$timeout){return(message("Failed to converge", group_2_text, " model prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))}


  both_groups_rhythmic <- ifelse(g1_model$rhythmic & g2_model$rhythmic, TRUE, FALSE)

  if(!both_groups_rhythmic){
    if(!g1_model$rhythmic & !g2_model$rhythmic){
      return(message("Both groups of data were arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }
    if(!g1_model$rhythmic){
      return(message(group_1_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }else{
      return(message(group_2_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }
  }
  n <- 0
  while(!comparison_model_success & !comparison_model_timeout){
    starting_params <- list(
      k=g1_model$k_estimate*2*stats::runif(1),
      k1=(g2_model$k_estimate - g1_model$k_estimate)*2*stats::runif(1),
      alpha=g1_model$alpha_estimate*2*stats::runif(1),
      alpha1=(g2_model$alpha_estimate - g1_model$alpha_estimate)*2*stats::runif(1),
      phi=g1_model$phi_estimate*2*stats::runif(1),
      phi1=random_start_phi1(p=(g2_model$phi_estimate - g1_model$phi_estimate))
    )
    fit.nls <- try({stats::nls(measure~k+k1*x_group+(alpha+alpha1*x_group)*cos(time_r-(phi+phi1*x_group)),
                               data = x,
                               start = starting_params,
                               stats::nls.control(maxiter = 100, minFactor = 1/10000)
                               )},
                   silent = TRUE)
    if (class(fit.nls) == "try-error") {
      n <- n + 1
      comparison_model_timeout <- ifelse(n>timeout_n, TRUE, FALSE)
    }else{
      nls_coefs <- extract_model_coefs(fit.nls)
      k_out <- stats::coef(fit.nls)[1]
      k1_out <- stats::coef(fit.nls)[2]
      k_out_p <- (summary(fit.nls)$coef)[1,4]
      k1_out_p <- (summary(fit.nls)$coef)[2,4]

      alpha_out <- stats::coef(fit.nls)[3]
      alpha1_out <- stats::coef(fit.nls)[4]
      alpha1_out_p <- (summary(fit.nls)$coef)[4,4]
      phi_out <- stats::coef(fit.nls)[5]
      phi1_out <- stats::coef(fit.nls)[6]
      phi1_out_p <- (summary(fit.nls)$coef)[6,4]

      comparison_model_success <- ifelse(nls_coefs['alpha','estimate']>0 & (nls_coefs['alpha','estimate'] + nls_coefs['alpha1','estimate']) > 0 & nls_coefs['phi1','estimate'] <pi & nls_coefs['phi1','estimate'] >-pi, 1, 0)
      comparison_model_timeout <- ifelse(n>timeout_n, TRUE, FALSE)
      n <- n + 1
    }
  }
  if(comparison_model_timeout){
    return(message("Both groups of data were rhythmic but the curve fitting procedure failed due to timing out. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
  }

  # Store the parameter estimate values in a named (numeric) vector for less verbose access
  V <- nls_coefs[, 'estimate']

  eq_1 <- function(time){V['k'] + V['alpha']*cos((2*pi/period)*time - V['phi'])}
  eq_2 <- function(time){V['k'] + V['k1'] + (V['alpha'] + V['alpha1'])*cos((2*pi/period)*time - (V['phi'] + V['phi1']))}
  fig_out <- ggplot2::ggplot(x, ggplot2::aes(time, measure)) +
    ggplot2::stat_function(ggplot2::aes(colour = group_1_text), fun = eq_1, size = 1) +
    ggplot2::stat_function(ggplot2::aes(colour = group_2_text), fun = eq_2, size = 1) +
    ggplot2::geom_point(ggplot2::aes(colour = group)) +
    ggplot2::scale_colour_manual(breaks = c(group_1_text, group_2_text),
                                 values = c("blue", "red")) +
    ggplot2::labs(colour = 'Legend',
                  x = "time (hours)")+
    ggplot2::xlim(min(floor(x$time/period) * period),
                  max(ceiling(x$time/period) * period))

  # Adjust phi_out so that -pi < phi_out < pi
  if(V['phi'] > pi){
    while(V['phi'] > pi){
      V['phi'] <- V['phi'] - 2*pi
    }
  }
  if(V['phi'] < -pi){
    while(V['phi'] < -pi){
      V['phi'] <- V['phi'] + 2*pi
    }
  }
  baseline_diff_abs <- V['k']
  baseline_diff_pc <- ((V['k'] + V['k1'])/V['k'])*100 - 100
  amplitude_diff_abs <- V['alpha1']
  amplitude_diff_pc <-  ((V['alpha']+V['alpha1'])/V['alpha'])*100 - 100
  g1_peak_time <- V['phi']*period/(2*pi)
  g2_peak_time <- (V['phi']+V['phi1'])*period/(2*pi)
  while(g1_peak_time > period | g1_peak_time < 0){
    if(g1_peak_time > period){
      g1_peak_time <- g1_peak_time - period
    }
    if(g1_peak_time<0){
      g1_peak_time <- g1_peak_time + period
    }
  }
  while(g2_peak_time>period| g2_peak_time < 0){
    if(g2_peak_time>period){
      g2_peak_time <- g2_peak_time - period
    }
    if(g2_peak_time<0){
      g2_peak_time <- g2_peak_time + period
    }
  }
  peak_time_diff <- V['phi1']*period/(2*pi)

  output_parms <- data.frame(parameter = c("Both groups were rhythmic",
                                           paste("Presence of rhythmicity (p-value) for ", group_1_text, sep = ""),
                                           paste("Presence of rhythmicity (p-value) for ", group_2_text, sep = ""),
                                           paste(group_1_text, " mesor estimate", sep = ""),
                                           paste(group_2_text, " mesor estimate", sep = ""),
                                           "Mesor difference estimate",
                                           "P-value for mesor difference",
                                           paste(group_1_text, " amplitude estimate", sep = ""),
                                           paste(group_2_text, " amplitude estimate", sep = ""),
                                           "Amplitude difference estimate",
                                           "P-value for amplitude difference",
                                           paste(group_1_text, " peak time", sep = ""),
                                           paste(group_2_text, " peak time", sep = ""),
                                           "Phase difference estimate",
                                           "P-value for difference in phase"),

                             value = c(both_groups_rhythmic, g1_model$alpha_p, g2_model$alpha_p, V['k'], (V['k'] + V['k1']), V['k1'],
                                       nls_coefs['k1', 'p_value'], V['alpha'], V['alpha'] + V['alpha1'], V['alpha1'], nls_coefs['alpha1', 'p_value'],
                                       g1_peak_time, g2_peak_time, peak_time_diff, nls_coefs['phi1', 'p_value']))

  return(list(plot=fig_out, table=output_parms, fit=fit.nls))
}
