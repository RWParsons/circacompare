make_data <- function(hours_diff){
  set.seed(123)
  g1 <- data.frame(time = c(),
                   measure = c())
  g2 <- data.frame(time = c(),
                   measure = c())
  for(i in 1:24){
    x1 <- (i*2*pi)/24
    x2 <- ((i+hours_diff)*2*pi)/24
    measure_1 <- sin(x1)*10
    measure_1 <- rnorm(1, measure_1, 4)
    measure_2 <- sin(x2)*14 + 3
    measure_2 <- rnorm(1, measure_2, 2)
    g1 <- rbind(g1, c(i,signif(measure_1,digits = 2)))
    g2 <- rbind(g2, c(i,signif(measure_2,digits = 2)))
  }
  colnames(g1) <- c("time", "measure")
  colnames(g2) <- c("time", "measure")
  g1$group <- "g1"
  g2$group <- "g2"
  df <- rbind(g1,g2)
}

circacompare <- function(x,
                          col_time,
                          col_group,
                          col_outcome,
                          period = 24,
                          alpha_threshold = 0.05,
                          timeout_n = 10000){
  
  if(!"ggplot2" %in% installed.packages()[, "Package"]){
    return(message("Please install 'ggplot2'"))
  }
  

  
  library(ggplot2)
  colnames(x)[agrep(col_group, colnames(x))] <- "group"

  if(length(levels(as.factor(x$group)))!=2){
    return(message("Your grouping variable had more or less than 2 levels! \nThis function is used to compare two groups of data. \nTo avoid me having to guess, please send data with only two possible values in your grouping variable to this function."))
  }
  
  group_1_text <- levels(as.factor(x$group))[1]
  group_2_text <- levels(as.factor(x$group))[2]
  colnames(x)[agrep(col_time, colnames(x))] <- "time"
  colnames(x)[agrep(col_outcome, colnames(x))] <- "measure"
  
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

  while(g1_success !=1){
    g1_alpha_start <- runif(1)*1000
    g1_phi_start <- runif(1)*6.15 - 3.15
    fit.nls_group_1 <- nls(measure~k + alpha*cos(time_r-phi),
                           data = dat_group_1,
                           start = list(k=1,alpha=g1_alpha_start,phi=g1_phi_start))
    g1_alpha_out <- summary(fit.nls_group_1)$coef[2,1]
    g1_alpha_p <- summary(fit.nls_group_1)$coef[2,4]
    g1_phi_out <- summary(fit.nls_group_1)$coef[3,1]
    g1_success <- ifelse(g1_alpha_out > 0,1,0)
  }
  while(g2_success !=1){
    g2_alpha_start <- runif(1)*1000
    g2_phi_start <- runif(1)*6.15 - 3.15
    fit.nls_group_2 <- nls(measure~k + alpha*cos(time_r-phi),
                           data = dat_group_2,
                           start = list(k=1,alpha=g2_alpha_start,phi=g2_phi_start))
    g2_alpha_out <- summary(fit.nls_group_2)$coef[2,1]
    g2_alpha_p <- summary(fit.nls_group_2)$coef[2,4]
    g2_phi_out <- summary(fit.nls_group_2)$coef[3,1]
    g2_success <- ifelse(g2_alpha_out > 0,1,0)
  }

  g1_rhythmic <- ifelse(g1_alpha_p < alpha_threshold, TRUE, FALSE)
  g2_rhythmic <- ifelse(g2_alpha_p < alpha_threshold, TRUE, FALSE)
  both_groups_rhythmic <- ifelse(g1_rhythmic ==TRUE & g2_rhythmic==TRUE, TRUE, FALSE)

  if(both_groups_rhythmic == TRUE){
    while(comparison_model_success == 0 & comparison_model_timeout == FALSE){
      alpha_in <- g1_alpha_out
      alpha1_in <- runif(1)*10 -5
      phi_in <- g1_phi_out
      phi1_in <- runif(1)*2*pi - pi
      fit.nls <- nls(measure~k+k1*x_group+(alpha+alpha1*x_group)*cos(time_r-(phi+phi1*x_group)),
                     data = x,
                     start = list(k=1, k1=0, alpha=alpha_in, alpha1=alpha1_in, phi=phi_in, phi1=phi1_in),
                     nls.control(maxiter = 100, minFactor = 1/10000, warnOnly = TRUE))

      k_out <- coef(fit.nls)[1]
      k1_out <- coef(fit.nls)[2]
      k_out_p <- (summary(fit.nls)$coef)[1,4]
      k1_out_p <- (summary(fit.nls)$coef)[2,4]

      alpha_out <- coef(fit.nls)[3]
      alpha1_out <- coef(fit.nls)[4]
      alpha1_out_p <- (summary(fit.nls)$coef)[4,4]

      phi_out <- coef(fit.nls)[5]
      phi1_out <- coef(fit.nls)[6]
      phi1_out_p <- (summary(fit.nls)$coef)[6,4]

      comparison_model_success <- ifelse(alpha_out>0 & (alpha_out + alpha1_out) > 0 & phi1_out <pi & phi1_out >-pi, 1, 0)
      n <- n + 1
      comparison_model_timeout <- ifelse(n>timeout_n, TRUE, FALSE)
    }
    
    if(comparison_model_timeout == TRUE){
      return(message("Both groups of data were rhythmic but the curve fitting procedure failed due to timing out. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
      }
    #loop curve fitting process (all data) until outputs are appropriate, or until looped more times than timeout_n
    if(comparison_model_timeout == FALSE){
      eq_1 <- function(time){k_out + alpha_out*cos((2*pi/period)*time - phi_out)}
      eq_2 <- function(time){k_out + k1_out + (alpha_out + alpha1_out)*cos((2*pi/period)*time - (phi_out + phi1_out))}
      
      fig_out <- ggplot2::ggplot(x, aes(time, measure)) +
        stat_function(fun = eq_1, colour = "deep sky blue", size=1) +
        stat_function(fun = eq_2, colour = "red", size=1) +
        geom_point(aes(colour = group)) + 
        scale_colour_manual(breaks = c(group_1_text, group_2_text),
                            values = c("deep sky blue", "red")) +
        xlab("time (hours)")
      
    #if the nls was successful, create a graph to plot the data as well as curves of best fit, 'fig_out'
  }
  if(both_groups_rhythmic==TRUE & comparison_model_success==1){
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
    g1_peak_time <- phi_out*24/(2*pi)
    g2_peak_time <- (phi_out+phi1_out)*24/(2*pi)
    while(g1_peak_time >24 | g1_peak_time < 0){
      if(g1_peak_time >24){
        g1_peak_time <- g1_peak_time - 24
      }
      if(g1_peak_time<0){
        g1_peak_time <- g1_peak_time + 24
      }
    }
    while(g2_peak_time >24| g2_peak_time <0){
      if(g2_peak_time>24){
        g2_peak_time <- g2_peak_time - 24
      }
      if(g2_peak_time<0){
        g2_peak_time <- g2_peak_time + 24
      }
    }
    peak_time_diff <- phi1_out*24/(2*pi)
  }

  if(comparison_model_timeout == TRUE | both_groups_rhythmic==FALSE){
    k_out <- NA
    k_out_p <- NA
    k1_out <- NA
    k1_out_p <- NA
    alpha_out <- NA
    alpha1_out <- NA
    alpha1_out_p <- NA
    phi_out <- NA
    phi1_out <- NA
    phi1_out_p <- NA
    baseline_diff_abs <- NA
    baseline_diff_pc <- NA
    amplitude_diff_abs <- NA
    amplitude_diff_pc <- NA
    g1_peak_time <- NA
    g2_peak_time <- NA
    peak_time_diff <- NA

    fig_out <- ggplot2::ggplot(x, aes(time, measure)) +
      geom_point(aes(colour = group))+
      scale_colour_manual(breaks = c(group_1_text, group_2_text),
                          values = c("deep sky blue", "red"))
  }
  output_parms <- data.frame(parameter = c("both_groups_rhythmic", 
                                           paste(group_1_text, "_rhythmicity_p", sep = ""), 
                                           paste(group_2_text, "_rhythmicity_p", sep = ""),
                                           paste(group_1_text, "_MESOR_estimate", sep = ""),
                                           "MESOR_difference_estimate", 
                                           "MESOR_difference_p", 
                                           paste(group_1_text, "_amplitude_estimate", sep = ""), 
                                           "amplitude_difference_estimate",
                                           "amplitude_difference_p", 
                                           paste(group_1_text, "_peak_time", sep = ""), 
                                           "phase_difference_estimate", 
                                           "phase_difference_p"),
                             value = c(both_groups_rhythmic, g1_alpha_p, g2_alpha_p, k_out, k1_out, 
                                       k1_out_p, alpha_out, alpha1_out, alpha1_out_p,
                                       g1_peak_time, peak_time_diff, phi1_out_p))


  if(exists("fig_out")){
    return(list(fig_out, output_parms))
  }
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
}
