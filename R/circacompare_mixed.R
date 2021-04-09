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
#' @param randomeffects which rhythmic parameters to allow random effects. The default is to include all rhythmic parameters.
#' @param period The period of the rhythm. For circadian rhythms, leave this as the default value, \code{24}.
#' @param alpha_threshold The level of alpha for which the presence of rhythmicity is considered. Default is to \code{0.05}.
#' @param nlme_control A list of control values for the estimation algorithm to replace the default values returned by the function nlme::nlmeControl. Defaults to an empty list.
#' @param nlme_method A character string. If "REML" the model is fit by maximizing the restricted log-likelihood. If "ML" the log-likelihood is maximized. Defaults to "REML".
#' @param verbose An optional logical value. If \code{TRUE} information on the evolution of the iterative algorithm is printed. Default is \code{FALSE}.
#' @param timeout_n The upper limit for the model fitting attempts. Default is \code{10000}.
#'
#' @return list
#' @export
#'
#' @examples
#' # Generate some data with within-id correlation for phase-shift (phi1)
#' mixed_data <- function(n){
#'     counter <- 1
#'     for(i in 1:n){
#'         x <- make_data(k1=0, alpha=0, phi1=rnorm(1, 6, 1))
#'         x$id <- counter
#'         counter <- counter + 1
#'         if(i==1){res <- x}else{res <- rbind(res, x)}
#'     }
#'     return(res)
#' }
#'
#' set.seed(42)
#' df <- mixed_data(n=50)
#'
#' # Use a mixed model with a random effect only on phase-relevant parameters as
#' # to not overparametrize the model.
#' out <- circacompare_mixed(x = df, col_time = "time", col_group = "group",
#'                           col_outcome = "measure", col_id = "id",
#'                           randomeffects = c("phi", "phi1"))
#' out
circacompare_mixed <- function(x,
                               col_time,
                               col_group,
                               col_outcome,
                               col_id,
                               randomeffects = c("k", "alpha", "phi"),
                               period = 24,
                               alpha_threshold = 0.05,
                               nlme_control = list(),
                               nlme_method = "REML",
                               verbose = FALSE,
                               timeout_n = 10000){


  if(!requireNamespace("ggplot2", quietly = TRUE)){
    return(message("Please install 'ggplot2'"))
  }

  if(!requireNamespace("nlme", quietly = TRUE)){
    return(message("Please install 'nlme'"))
  }

  x <- x[c(col_time, col_group, col_id, col_outcome)]
  colnames(x) <- c("time", "group", "id", "measure")

  if(length(setdiff(randomeffects, c("k", "k1", "alpha", "alpha1", "phi", "phi1"))) != 0){
    return(message('"randomeffects" should only include the names of parameters\nthat represent rhythmic characteristics in the model.\nThey should be a subset of, or equal to c("k", "k1", "alpha", "alpha1", "phi", "phi1")'))
  }

  if(length(randomeffects) == 0){
    return(message("If you do not want to include any random effects, than you ought to use 'circacompare' rather than 'circacompare_mixed'"))
  }

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

  dat_group_1 <- x[x$group == group_1_text,]
  dat_group_2 <- x[x$group == group_2_text,]

  g1_model <- model_each_group(data=dat_group_1, type="nlme",
                               args=list(
                                 timeout_n=timeout_n,
                                 alpha_threshold=alpha_threshold,
                                 randomeffects=randomeffects,
                                 nlme_method=nlme_method,
                                 nlme_control=nlme_control,
                                 verbose=verbose
                               ))
  if(g1_model$timeout){return(message("Failed to converge", group_1_text, " model prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))}

  g2_model <- model_each_group(data=dat_group_2, type="nlme",
                               args=list(
                                 timeout_n=timeout_n,
                                 alpha_threshold=alpha_threshold,
                                 randomeffects=randomeffects,
                                 nlme_method=nlme_method,
                                 nlme_control=nlme_control,
                                 verbose=verbose
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
  comparison_model_success <- FALSE
  comparison_model_timeout <- FALSE
  while(!comparison_model_success & !comparison_model_timeout){
    starting_params <- c(
      k=g1_model$k_estimate*2*stats::runif(1),
      k1=(g2_model$k_estimate - g1_model$k_estimate)*2*stats::runif(1),
      alpha=g1_model$alpha_estimate*2*stats::runif(1),
      alpha1=(g2_model$alpha_estimate - g1_model$alpha_estimate)*2*stats::runif(1),
      phi=g1_model$phi_estimate*2*stats::runif(1),
      phi1=random_start_phi1(p=(g2_model$phi_estimate - g1_model$phi_estimate))
    )
    randomeffects_formula <- stats::formula(paste(paste0(randomeffects, collapse="+"), "~ 1 | id"))
    fit.nlme <- try({nlme::nlme(measure~k+k1*x_group+(alpha+alpha1*x_group)*cos(time_r-(phi+phi1*x_group)),
                                random = randomeffects_formula,
                                fixed = k+k1+alpha+alpha1+phi+phi1~1,
                                data = x,
                                start = starting_params,
                                method=nlme_method,
                                control = nlme_control,
                                verbose = verbose)},
                    silent = ifelse(verbose, FALSE, TRUE))

    if ("try-error" %in% class(fit.nlme)) {
      n <- n + 1
    }
    else{
      nlme_coefs <- extract_model_coefs(fit.nlme)
      comparison_model_success <- ifelse(nlme_coefs['alpha', 'estimate']>0 & (nlme_coefs['alpha', 'estimate'] + nlme_coefs['alpha1', 'estimate']) > 0 & nlme_coefs['phi1', 'estimate'] < pi & nlme_coefs['phi1', 'estimate'] >- pi, TRUE, FALSE)
      comparison_model_timeout <- ifelse(n>timeout_n, TRUE, FALSE)
      n <- n + 1
    }
  }
  if(comparison_model_timeout){
    return(message("Both groups of data were rhythmic but the curve fitting procedure failed due to timing out. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
  }

  # Store the parameter estimate values in a named (numeric) vector for less verbose access
  V <- nlme_coefs[, 'estimate']

  # Graph the model and data
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

  #  Adjust phi_out so that -pi < phi_out < pi
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
                                       nlme_coefs['k1', 'p_value'], V['alpha'], V['alpha'] + V['alpha1'], V['alpha1'], nlme_coefs['alpha1', 'p_value'],
                                       g1_peak_time, g2_peak_time, peak_time_diff, nlme_coefs['phi1', 'p_value']))

  return(list(plot=fig_out, table=output_parms, fit=fit.nlme))
}

