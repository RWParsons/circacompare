\name{circa_single}
\alias{circa_single}
\title{Analysis of circadian data}
\usage{
circa_single(x, col_time, col_outcome, period = 24, alpha_threshold = 0.05, timeout_n = 10000)
}
\arguments{
\item{x}{A data.frame containing circadian data.}

\item{col_time}{Character string which corresponds to the name of the variable in the dataset containing the time (numeric and in hours)}

\item{col_outcome}{Character string which corresponds to the name of the variable in the dataset containing the outcome.}

\item{period}{The period of the data (in hours).  For circadian data, use \code{period = 24}.}

\item{alpha_threshold}{This determines the alpha value for classifying presence of rhythmicity. For example, if \code{alpha_threshold = 0.05}, then the data will be classified as rhythmic if p < 0.05 for amplitude being greater than zero.}

\item{timeout_n}{A value which controls the upper limit for the nls fitting attempts.}
}
\description{
Takes circadian data in a single dataset and performs a nonlinear regression.  Values are estimated and given for mesor, amplitude and phase. Acommpanying the summary results, a ggplot and nls (nonlinear least squares regression) object is also returned.  The nls object can be used to extract p-values for the parameters difference from zero as well as the standard error for the given estimate. The data which is passed to \code{circa_single} does not necessarily have to include times only within hours 0 to 24. 
}
\examples{
# create an example data frame with a phase shift between groups of 6 hours.
df <- make_data(hours_diff = 6)
df <- df[df$group == "g1",]

# call circacompare.
out <- circacompare(x = df, col_time = "time", col_outcome = "measure")

# view the graph.
out[[1]]

# view the results.
out[[2]]
  
# view the nls object summary and confidence intervals for parameters.  To understand the parameters in this object, it would be useful to read the circacompare publication (https://doi.org/10.1093/bioinformatics/btz730)
summary(out[[3]])
nlstools::confint2(out[[3]])

}
