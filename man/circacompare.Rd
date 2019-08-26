\name{circacompare}
\alias{circacompare}
\title{'statistically compares two groups of circadian data'}
\usage{
circacompare(x, col_time, col_group, col_outcome, period = 24, alpha_threshold = 0.05, timeout_n = 10000)
}
\description{
Takes circadian data of two groups in a single dataset and tests each group for rhythmicity. If both groups are rhythmic, the two are compared. Statistical and estimated differences are given for mesor, amplitude and phase. Acommpanying the summary results, a ggplot object is also returned.
}
\examples{
df <- make_data(hours_diff = 6)
circacompare(x = df, col_time = "time", col_group = "group", col_outcome = "measure")
}