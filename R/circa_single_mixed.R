#' @title circa_single_mixed
#' @name circa_single_mixed
#'
#' @description \code{circa_single_mixed} is similar to \code{circa_single} but allows for some simple, user-specified random-effects on the rhythmic parameters of choice.
#'
#' @param x data.frame.  This is the data.frame which contains the rhythmic data in a tidy format.
#' @param col_time The name of the column within the data.frame, x, which contains time in hours at which the data were collected.
#' @param col_outcome The name of the column within the data.frame, x, which contains outcome measure of interest.
#' @param col_id The name of the column within the data.frame, \code{x}, which contains the identifying values for the random effect, such as \code{subject_id}.
#' @param randomeffects which rhythmic parameters to allow random effects. The default is to include all rhythmic parameters.
#' @param period The period of the rhythm. For circadian rhythms, leave this as the default value, \code{24}.
#' @param alpha_threshold The level of alpha for which the presence of rhythmicity is considered. Default is to \code{0.05}.
#' @param nlme_control A list of control values for the estimation algorithm to replace the default values returned by the function nlme::nlmeControl. Defaults to an empty list.
#' @param nlme_method A character string. If "REML" the model is fit by maximizing the restricted log-likelihood. If "ML" the log-likelihood is maximized. Defaults to "ML".
#' @param verbose An optional logical value. If \code{TRUE} information on the evolution of the iterative algorithm is printed. Default is \code{FALSE}.
#' @param timeout_n The upper limit for the model fitting attempts. Default is \code{10000}.
#' @param return_figure Whether or not to return a ggplot graph of the rhythm and cosine model.
#'
#' @return list
#' @export
#'
#' @examples
#' set.seed(42)
#' mixed_data <- function(n){
#'   counter <- 1
#'   for(i in 1:n){
#'       x <- make_data(k1=rnorm(1, 10, 2), alpha1=0, phi1=0)
#'       x$id <- counter
#'       counter <- counter + 1
#'       if(i==1){res <- x}else{res <- rbind(res, x)}
#'   }
#'   return(res)
#' }
#' df <- mixed_data(n=50)
#' out <- circa_single_mixed(x = df, col_time = "time", col_outcome = "measure",
#'                          col_id = "id", randomeffects = c("k"))
circa_single_mixed <- function (x,
                                col_time,
                                col_outcome,
                                col_id,
                                randomeffects = c("k", "alpha", "phi"),
                                period = 24,
                                alpha_threshold = 0.05,
                                nlme_control = list(),
                                nlme_method = "ML",
                                verbose = FALSE,
                                timeout_n = 10000,
                                return_figure = TRUE){
  if(!requireNamespace("ggplot2", quietly = TRUE) & return_figure){
    return(message("Please install 'ggplot2'"))
  }

  x <- x[c(col_time, col_outcome, col_id)]
  colnames(x) <- c("time", "measure", "id")

  if(length(setdiff(randomeffects, c("k", "alpha", "phi"))) != 0){
    return(message('"randomeffects" should only include the names of parameters\nthat represent rhythmic characteristics in the model.\nThey should be a subset of, or equal to c("k", "alpha", "phi")'))
  }

  if(length(randomeffects) == 0){
    return(message("If you do not want to include any random effects, than you ought to use 'circa_single' rather than 'circa_single_mixed'"))
  }

  if (!class(x$time) %in% c("numeric", "integer")) {
    return(message(paste("The time variable which you gave was a '",
                         class(x$time), "' \nThis function expects time to be given as hours and be of class 'integer' or 'numeric'.",
                         "\nPlease convert the time variable in your dataframe to be of one of these classes",
                         sep = "")))
  }

  if (!class(x$measure) %in% c("numeric", "integer")) {
    return(message(paste("The measure variable which you gave was a '",
                         class(x$measure), "' \nThis function expects measure to be number and be of class 'integer' or 'numeric'.",
                         "\nPlease convert the measure variable in your dataframe to be of one of these classes",
                         sep = "")))
  }
  x$time_r <- (x$time/24) * 2 * pi * (24/period)
  comparison_model_success <- 0
  comparison_model_timeout <- FALSE
  success <- 0
  n <- 0
  while (success != 1) {
    alpha_start <- (max(x$measure, na.rm = TRUE) - min(x$measure, na.rm = TRUE)) * stats::runif(1)
    phi_start <- stats::runif(1) * 6.15 - 3.15
    k_start <- mean(x$measure, na.rm = TRUE) * 2 * stats::runif(1)

    randomeffects_formula <- stats::formula(paste(paste0(randomeffects, collapse="+"), "~ 1 | id"))

    fit.nlme <- try({nlme::nlme(measure~k+alpha*cos(time_r-phi),
                                random = randomeffects_formula,
                                fixed = k+alpha+phi~1,
                                data = x,
                                start = c(k=k_start, alpha=alpha_start, phi=phi_start),
                                control = nlme_control,
                                method = nlme_method,
                                verbose = verbose)},
                    silent = ifelse(verbose, FALSE, TRUE)
    )
    if ("try-error" %in% class(fit.nlme)) {
      n <- n + 1
    }
    else {
      k_out <- summary(fit.nlme)$tTable[1, 1]
      alpha_out <- summary(fit.nlme)$tTable[2, 1]
      alpha_p <- summary(fit.nlme)$tTable[2, 5]
      phi_out <- summary(fit.nlme)$tTable[3, 1]
      success <- ifelse(alpha_out > 0 & phi_out >= 0 & phi_out <= 2 * pi, 1, 0)
      n <- n + 1
    }
    if (n >= timeout_n) {
      return(message("Failed to converge data prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
  }
  data_rhythmic <- ifelse(alpha_p < alpha_threshold, TRUE, FALSE)
  eq <- function(time) {
    k_out + alpha_out * cos((2 * pi/period) * time - phi_out)
  }

  if(return_figure == TRUE){
    if(data_rhythmic == TRUE) {
      fig_out <- ggplot2::ggplot(x, ggplot2::aes(time, measure)) +
        ggplot2::stat_function(fun = eq, size = 1) +
        ggplot2::geom_point() +
        ggplot2::xlim(min(floor(x$time/period) * period),
                      max(ceiling(x$time/period) * period)) +
        ggplot2::labs(subtitle = "Data is rhythmic", x = "time (hours)")
    }
    else{
      fig_out <- ggplot2::ggplot(x, ggplot2::aes(time, measure)) +
        ggplot2::geom_point() +
        ggplot2::xlim(min(floor(x$time/period) * period),
                      max(ceiling(x$time/period) * period)) +
        ggplot2::labs(subtitle = "Data is arrhythmic", x = "time (hours)")
    }
  }

  peak_time <- phi_out * period/(2 * pi)
  output_parms <- data.frame(mesor = k_out, amplitude = alpha_out,
                             amplitude_p = alpha_p, phase_radians = phi_out, peak_time_hours = phi_out *
                               period/(2 * pi))
  if(return_figure == TRUE){
    return(list(fit.nlme, output_parms, fig_out))
  }else{
    return(list(fit.nlme, output_parms))
  }

}
