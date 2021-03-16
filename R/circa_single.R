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
#'
#' @return list
#' @export
#'
#' @examples
#' df <- make_data(phi1 = 0)
#' circa_single(x = df, col_time = "time", col_outcome="measure")
circa_single <- function (x,
                          col_time,
                          col_outcome,
                          period = 24,
                          alpha_threshold = 0.05,
                          timeout_n = 10000,
                          return_figure = TRUE){
  if(!requireNamespace("ggplot2", quietly = TRUE) & return_figure){
    return(message("Please install 'ggplot2'"))
  }

  x <- x[c(col_time, col_outcome)]
  colnames(x) <- c("time", "measure")

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
    fit.nls <- try({
      stats::nls(measure ~ k + alpha * cos(time_r - phi),
                 data = x,
                 start = list(k = k_start, alpha = alpha_start, phi = phi_start))
    }, silent = TRUE)
    if (class(fit.nls) == "try-error") {
      n <- n + 1
    }
    else {
      k_out <- summary(fit.nls)$coef[1, 1]
      alpha_out <- summary(fit.nls)$coef[2, 1]
      alpha_p <- summary(fit.nls)$coef[2, 4]
      phi_out <- summary(fit.nls)$coef[3, 1]
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

  if(return_figure){
    if(data_rhythmic) {
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
  if(return_figure){
    return(list(fit.nls, output_parms, fig_out))
  }else{
    return(list(fit.nls, output_parms))
  }

}
