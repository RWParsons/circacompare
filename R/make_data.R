#' @title make_data
#' @name make_data
#'
#' @description Generate example circadian data with specified phase shift between groups
#'
#' @param k1 mesor difference between groups.
#' @param alpha1 amplitude difference between groups.
#' @param phi1 Amount of phase shift, in hours, between the two groups within the generated data
#' @param seed random seed for generating data.
#' @return data.frame
#' @export
#'
#' @examples
#' data <- make_data(k1=3, alpha1=4, phi1 = 6)
make_data <- function(k1=3, alpha1=4, phi1 = 6, seed=NULL){
  if(!is.null(seed)){set.seed(seed)}
  g1 <- data.frame(time = c(),
                   measure = c())
  g2 <- data.frame(time = c(),
                   measure = c())
  for(i in 1:24){
    x1 <- (i*2*pi)/24
    x2 <- ((i+phi1)*2*pi)/24
    measure_1 <- sin(x1)*10
    measure_1 <- stats::rnorm(1, measure_1, 2)
    measure_2 <- sin(x2)*(10 + alpha1) + k1
    measure_2 <- stats::rnorm(1, measure_2, 2)
    g1 <- rbind(g1, c(i,signif(measure_1,digits = 2)))
    g2 <- rbind(g2, c(i,signif(measure_2,digits = 2)))
  }
  colnames(g1) <- c("time", "measure")
  colnames(g2) <- c("time", "measure")
  g1$group <- "g1"
  g2$group <- "g2"
  df <- rbind(g1,g2)
}
