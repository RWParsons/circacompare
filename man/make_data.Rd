\name{make_data}
\alias{make_data}
\title{'makes example circadian data'}
\usage{
make_data(hours_diff)
}
\description{
Creates example dataset two groups of circadian data.  The phase difference between the two groups in the created data is specified by the user using the `hours_diff` argument.
}
\examples{
df <- make_data(6)
}