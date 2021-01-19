\name{circacompare_mixed}
\alias{circacompare_mixed}
\title{Statistically compares two groups of circadian data with a hierarchical structure}
\usage{
circacompare(x, col_time, col_group, col_outcome, col_id, period = 24,
			 alpha_threshold = 0.05, nlme_control = list(), verbose = FALSE, 
			 timeout_n = 10000))
}
\arguments{
\item{x}{A data.frame containing circadian data of two groups.}

\item{col_time}{Character string which corresponds to the name of the variable in the dataset containing the time (numeric and in hours)}

\item{col_group}{Character string which corresponds to the name of the variable in the dataset used for grouping the data.  There should only be two levels for this variable.}

\item{col_outcome}{Character string which corresponds to the name of the variable in the dataset containing the outcome.}

\item{col_id}{Character string which corresponds to the name of the variable in the dataset on which the hierarchical structure is placed, such as the subject ID}

\item{period}{The period of the data (in hours).  For circadian data, use \code{period = 24}.}

\item{alpha_threshold}{This determines the alpha value for classifying presence of rhythmicity. For example, if \code{alpha_threshold = 0.05}, then the data will be classified as rhythmic if p < 0.05 for amplitude being greater than zero.}

\item{nlme_control}{This is a \code{nlmeControl} object which gets passed to the nonlinear mixed model fitting function, \code{nlme}. See \code{?nlme::nlmeControl}}

\item{timeout_n}{A value which controls the upper limit for the model fitting attempts.}
}
\description{
Takes circadian data of two groups in a single dataset and tests each group for rhythmicity. If both groups are rhythmic, the two are compared.  A heirarchical structure is applied byusing a mixed model.  The mixed model allows for a \code{col_id} level random effect on each of the rhythmic parameters: mesor, amplitude and phase. Statistical and estimated differences (fixed effects) are given for mesor, amplitude and phase. Acommpanying the summary results, a ggplot and nlme (nonlinear mixed-effects model) object is also returned.  The data which is passed to \code{circacompare_mixed} does not necessarily have to include times only within hours 0 to 24. 
}
\examples{
# create an example data frame with a phase shift between groups of 6 hours.
mixed_data <- function(n){
  counter <- 1
  for(i in 1:n){
    x <- make_data(rnorm(1, 6, 1))
    x$id[x$group=="g1"] <- counter
    x$id[x$group=="g2"] <- counter 
    counter <- counter + 1
    if(i==1){
      res <- x
    }else{
      res <- rbind(res, x)
    }
  }
  return(res)
}
df <- mixed_data(50)

# call circacompare.
out <- circacompare_mixed(x=df, col_time="time", col_group="group", col_outcome="measure", col_id="id")

# view the graph.
out[[1]]

# view the results.
out[[2]]
  
# view the nlme object summary and confidence intervals for the fixed effects.  To understand the parameters in this object, it may be useful to read the circacompare publication (https://doi.org/10.1093/bioinformatics/btz730).
summary(out[[3]])
intervals(out[[3]], which="fixed")


}
