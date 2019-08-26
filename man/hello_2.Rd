\name{circacompare}
\alias{circacompare}
\title{Illustration of crayon colors}
\usage{
circacompare(method2order = c("hsv", "cluster"), cex = 0.6,
             mar = rep(0.1, 4))
}
\arguments{
\item{method2order}{method to order colors (\code{"hsv"} or \code{"cluster"})}

\item{cex}{character expansion for the text}

\item{mar}{margin parameters; vector of length 4 (see \code{\link[graphics]{par}})}
}
\value{
None
}
\description{
Creates a plot of the crayon colors in \code{\link{brocolors}}
}
\examples{
circacompare()
}
