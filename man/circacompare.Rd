\name{circacompare}
\alias{circacompare}
\title{Statistically compares two groups of circadian data}
\usage{
circacompare(x, col_time, col_group, col_outcome,
             period = 24, alpha_threshold = 0.05, timeout_n = 10000)
}
\arguments{
\item{x}{A data.frame containing circadian data of two groups.}

\item{col_time}{Character string which corresponds to the name of the variable in the dataset containing the time (numeric and in hours)}

\item{col_group}{Character string which corresponds to the name of the variable in the dataset used for grouping the data.  There should only be two levels for this variable.}

\item{col_outcome}{Character string which corresponds to the name of the variable in the dataset containing the outcome.}

\item{period}{The period of the data (in hours).  For circadian data, use \code{period = 24}.}

\item{alpha_threshold}{This determines the alpha value for classifying presence of rhythmicity. For example, if \code{alpha_threshold = 0.05}, then the data will be classified as rhythmic if p < 0.05 for amplitude being greater than zero.}

\item{timeout_n}{A value which controls the upper limit for the nls fitting attempts.}
}
\description{
Takes circadian data of two groups in a single dataset and tests each group for rhythmicity. If both groups are rhythmic, the two are compared. Statistical and estimated differences are given for mesor, amplitude and phase. Acommpanying the summary results, a ggplot object is also returned.
}
\examples{
# create an example data frame with a phase shift between groups of 6 hours.
df <- make_data(hours_diff = 6)

# call circacompare.
out <- circacompare(x = df, col_time = "time", col_group = "group", col_outcome = "measure")

# View the graph.
out[[1]]

# View the results.
out[[2]]
}
