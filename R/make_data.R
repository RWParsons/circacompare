#' @title make_data
#' @name make_data
#'
#' @description Generate example circadian data with specified phase shift between groups
#'
#' @param k mesor of group 1.
#' @param k1 change in mesor in group 2 from group 1.
#' @param alpha amplitude rhythm for group 1.
#' @param alpha1 change in amplitude in group 2 from group 1
#' @param phi phase of rhythm, in radian-hours, in group 1.
#' @param phi1 change in phase, in radian-hours, in group 2 from group 1
#' @param hours the number of hours/datapoints to sample.
#' @param noise_sd the standard deviation of the noise term.
#' @param seed random seed for generating data.
#' @return data.frame
#' @export
#'
#' @examples
#' data <- make_data(k1=3, alpha1=4, phi1 = 6)
make_data <- function(k=0, k1=3, alpha=10, alpha1=4, phi=0, phi1 = 3.15, hours=48, noise_sd=0.1, seed=NULL){
  if(!is.null(seed)){set.seed(seed)}
  g1 <- data.frame(time = rep(NA, hours),
                   measure = rep(NA, hours))
  g2 <- data.frame(time = rep(NA, hours),
                   measure = rep(NA, hours))
  V <- c(k=k, k1=k1, alpha=alpha, alpha1=alpha1, phi=phi, phi1=phi1)
  eq_expression <- create_formula(grouped_params = c("k", "alpha", "phi"))$f_equation
  eval(parse(text=eq_expression$g1))
  eval(parse(text=eq_expression$g2))
  g1$time <- 1:hours
  g1$measure <- eq_1(g1$time) + stats::rnorm(n=hours, mean=0, sd=noise_sd)
  g1$group <- "g1"
  g2$time <- 1:hours
  g2$measure <- eq_2(g2$time) + stats::rnorm(n=hours, mean=0, sd=noise_sd)
  g2$group <- "g2"

  return(rbind(g1, g2))
}
