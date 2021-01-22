#' circacompare_mixed
#' @description \code{circacompare_mixed} is similar to \code{circacompare} but allows for some simple, user-specified random-effects on the rhythmic parameters of choice.
#'
#' @param x \code{data.frame}.  This is the data.frame which contains the rhythmic data for two groups in a tidy format.
#' @param col_time The name of the column within the data.frame, \code{x}, which contains time in hours at which the data were collected.
#' @param col_group The name of the column within the data.frame, \code{x}, which contains the grouping variable.  This should only have two levels.
#' @param col_outcome The name of the column within the data.frame, \code{x}, which contains outcome measure of interest.
#' @param col_id The name of the column within the data.frame, \code{x}, which contains the identifying values for the random effect, such as \code{subject_id}.
#' @param randomeffects which rhythmic parameters to allow random effects. The default is to include all rhythmic parameters.
#' @param period The period of the rhythm. For circadian rhythms, leave this as the default value, \code{24}.
#' @param alpha_threshold The level of alpha for which the presence of rhythmicity is considered. Default is to {0.05}.
#' @param nlme_control The list of control variables which is passed to the mixed model fitting function.  Create this using \code{nlme::nlmeControl(...)}.  This is left empty by default but you may need to adjust this if you're having convergence problems.
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
                               randomeffects = c("k", "k1", "alpha", "alpha1", "phi", "phi1"),
                               period = 24,
                               alpha_threshold = 0.05,
                               nlme_control = list(),
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

  comparison_model_success <- FALSE
  comparison_model_timeout <- FALSE
  g1_success <- FALSE
  g2_success <- FALSE
  g1_alpha_p <- NA
  g2_alpha_p <- NA
  n <- 0
  dat_group_1 <- x[x$group == group_1_text,]
  dat_group_2 <- x[x$group == group_2_text,]

  while(!g1_success){
    g1_alpha_start <- (max(dat_group_1$measure, na.rm = TRUE) - min(dat_group_1$measure, na.rm = TRUE)) * stats::runif(1)
    g1_phi_start <- stats::runif(1)*6.15 - 3.15
    g1_k_start <- mean(dat_group_1$measure, na.rm = TRUE)*2*stats::runif(1)

    fit.nlme_group_1 <- try({stats::nls(measure~k + alpha*cos(time_r-phi),
                                 data = dat_group_1,
                                 start = list(k=g1_k_start,alpha=g1_alpha_start,phi=g1_phi_start))},
                            silent = TRUE)
    if(class(fit.nlme_group_1) == "try-error"){
      n <- n + 1
    }
    else{
      g1_k_out <- summary(fit.nlme_group_1)$coef[1,1]
      g1_alpha_out <- summary(fit.nlme_group_1)$coef[2,1]
      g1_alpha_p <- summary(fit.nlme_group_1)$coef[2,4]
      g1_phi_out <- summary(fit.nlme_group_1)$coef[3,1]
      g1_success <- ifelse(g1_alpha_out > 0 & g1_phi_out > 0 & g1_phi_out < 2*pi ,TRUE,FALSE)
      n <- n + 1
      # if(g1_success){
      #   message("g1 success")
      #   return(fit.nlme_group_1)
      # }

    }
    if(n >= timeout_n){
      return(message("Failed to converge group 1 data prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
  }
  n <- 0
  while(!g2_success){
    g2_alpha_start <- (max(dat_group_2$measure, na.rm = TRUE) - min(dat_group_2$measure, na.rm = TRUE)) * stats::runif(1)
    g2_phi_start <- stats::runif(1)*6.15 - 3.15
    g2_k_start <- mean(dat_group_2$measure, na.rm = TRUE)*2*stats::runif(1)

    fit.nlme_group_2 <- try({stats::nls(measure~k + alpha*cos(time_r-phi),
                                 data = dat_group_2,
                                 start = list(k=g2_k_start,alpha=g2_alpha_start,phi=g2_phi_start))},
                            silent = TRUE)
    if(class(fit.nlme_group_2) == "try-error"){
      n <- n + 1
    }
    else{
      g2_k_out <- summary(fit.nlme_group_2)$coef[1,1]
      g2_alpha_out <- summary(fit.nlme_group_2)$coef[2,1]
      g2_alpha_p <- summary(fit.nlme_group_2)$coef[2,4]
      g2_phi_out <- summary(fit.nlme_group_2)$coef[3,1]
      g2_success <- ifelse(g2_alpha_out > 0 & g2_phi_out > 0 & g2_phi_out < 2*pi, TRUE, FALSE)
      n <- n + 1
      # if(g2_success){
      #   message("g2 success")
      #   return(fit.nlme_group_2)
      # }
    }
    if(n >= timeout_n){
      return(message("Failed to converge group 2 data prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
  }

  g1_rhythmic <- ifelse(g1_alpha_p < alpha_threshold, TRUE, FALSE)
  g2_rhythmic <- ifelse(g2_alpha_p < alpha_threshold, TRUE, FALSE)
  both_groups_rhythmic <- ifelse(g1_rhythmic ==TRUE & g2_rhythmic==TRUE, TRUE, FALSE)

  if(both_groups_rhythmic==FALSE){
    if(g1_rhythmic == FALSE & g2_rhythmic == FALSE){
      return(message("Both groups of data were arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }
    if(g1_rhythmic == FALSE){
      return(message(group_1_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }else{
      return(message(group_2_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }
  }
  n <- 0
  if(both_groups_rhythmic){
    while(!comparison_model_success & !comparison_model_timeout){
      alpha_in <- g1_alpha_out*2*stats::runif(1)
      alpha1_in <- (g2_alpha_out - g1_alpha_out)*2*stats::runif(1)
      phi_in <- g1_phi_out*2*stats::runif(1)

      phi1_in <- g2_phi_out - g1_phi_out
      phi1_in <- stats::rnorm(1, phi1_in, 2)
      if(abs(phi1_in) > 2*pi){
        if(phi1_in < 0){
          phi1_in <- phi1_in + 2*pi
        }else{
          phi1_in <- phi1_in - 2*pi
        }
      }


      k_in <- g1_k_out*2*stats::runif(1)
      k1_in <- (g2_k_out - g1_k_out)*2*stats::runif(1)

      randomeffects_formula <- stats::formula(paste(paste0(randomeffects, collapse="+"), "~ 1 | id"))

      fit.nlme <- try({nlme::nlme(measure~k+k1*x_group+(alpha+alpha1*x_group)*cos(time_r-(phi+phi1*x_group)),
                                  random = randomeffects_formula,
                                  fixed = k+k1+alpha+alpha1+phi+phi1~1,
                                  data = x,
                                  start = c(k=k_in, k1=k1_in, alpha=alpha_in, alpha1=alpha1_in, phi=phi_in, phi1=phi1_in),
                                  # control = nlme_control,
                                  verbose = verbose)},
                      silent = ifelse(verbose, FALSE, TRUE)
      )

      if ("try-error" %in% class(fit.nlme)) {
        n <- n + 1
      }
      else{
        k_out <- summary(fit.nlme)$tTable[1,1]
        k1_out <- summary(fit.nlme)$tTable[2,1]
        k_out_p <- summary(fit.nlme)$tTable[1,5]
        k1_out_p <- summary(fit.nlme)$tTable[2,5]

        alpha_out <- summary(fit.nlme)$tTable[3,1]
        alpha1_out <- summary(fit.nlme)$tTable[4,1]
        alpha1_out_p <- summary(fit.nlme)$tTable[4,5]

        phi_out <- summary(fit.nlme)$tTable[5,1]
        phi1_out <- summary(fit.nlme)$tTable[6,1]
        phi1_out_p <- summary(fit.nlme)$tTable[6,5]

        comparison_model_success <- ifelse(alpha_out>0 & (alpha_out + alpha1_out) > 0 & phi1_out < pi & phi1_out >- pi, TRUE, FALSE)
        comparison_model_timeout <- ifelse(n>timeout_n, TRUE, FALSE)
        n <- n + 1
      }
    }

    if(comparison_model_timeout){
      return(message("Both groups of data were rhythmic but the curve fitting procedure failed due to timing out. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
    #loop curve fitting process (all data) until outputs are appropriate, or until looped more times than timeout_n
    if(!comparison_model_timeout){
      eq_1 <- function(time){k_out + alpha_out*cos((2*pi/period)*time - phi_out)}
      eq_2 <- function(time){k_out + k1_out + (alpha_out + alpha1_out)*cos((2*pi/period)*time - (phi_out + phi1_out))}

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

    }
  }
  if(both_groups_rhythmic& comparison_model_success){
    if(phi_out > pi){
      while(phi_out > pi){
        phi_out <- phi_out - 2*pi
      }
    }
    if(phi_out < -pi){
      while(phi_out < -pi){
        phi_out <- phi_out + 2*pi
      }
    }#adjust phi_out so that -pi < phi_out < pi
    baseline_diff_abs <- k1_out
    baseline_diff_pc <- ((k_out + k1_out)/k_out)*100 - 100
    amplitude_diff_abs <- alpha1_out
    amplitude_diff_pc <-  ((alpha_out+alpha1_out)/alpha_out)*100 - 100
    g1_peak_time <- phi_out*period/(2*pi)
    g2_peak_time <- (phi_out+phi1_out)*period/(2*pi)
    while(g1_peak_time > period | g1_peak_time < 0){
      if(g1_peak_time > period){
        g1_peak_time <- g1_peak_time - period
      }
      if(g1_peak_time<0){
        g1_peak_time <- g1_peak_time + period
      }
    }
    while(g2_peak_time>period| g2_peak_time <0){
      if(g2_peak_time>period){
        g2_peak_time <- g2_peak_time - period
      }
      if(g2_peak_time<0){
        g2_peak_time <- g2_peak_time + period
      }
    }
    peak_time_diff <- phi1_out*period/(2*pi)
  }
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
                             value = c(both_groups_rhythmic, g1_alpha_p, g2_alpha_p, k_out, (k_out + k1_out), k1_out,
                                       k1_out_p, alpha_out, alpha_out + alpha1_out, alpha1_out, alpha1_out_p,
                                       g1_peak_time, g2_peak_time, peak_time_diff, phi1_out_p))


  if(exists("fig_out")){
    return(list(fig_out, output_parms, fit.nlme))
  }
}
