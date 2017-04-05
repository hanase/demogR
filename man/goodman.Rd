\name{goodman}
\alias{goodman}
\docType{data}
\title{ Demographic data from Venezuela, Madagascar and the United
States in the late 1960s}
\description{
  Aggregate vital registration data from Venezuela (1965), Madagascar
  (1966), and the United States (1967) used in Goodman, Keyfitz, and
  Pullam (1974).
}
\usage{data(goodman)}
\format{
  A data frame with 19 observations on the following 10 variables.
  \describe{
    \item{\code{age}}{ age classes }
    \item{\code{ven.nKx}}{ mid-year population structure for Venezuela}
    \item{\code{ven.nDx}}{ enumerated deaths for Venezuela}
    \item{\code{ven.bx}}{ enumerated births for Venezuela}
    \item{\code{mad.nKx}}{mid-year population structure for Madagascar}
    \item{\code{mad.nDx}}{ enumerated deaths for Madagascar }
    \item{\code{mad.bx}}{ enumerated births for Madagascar}
    \item{\code{usa.nKx}}{ mid-year population structure for the
    United States}
    \item{\code{usa.nDx}}{ enumerated deaths for the United States}
    \item{\code{usa.bx}}{ enumerated births for the United States}
  }
}
\details{
  Vital event data tabulated by Keyfitz and Flieger (1990) and used by
  Goodman et al. (1974) to illustrate their kinship frequency
  calculations.  The three populations represent tremendous diversity
  in age-specific demographic schedules.  Venezuela (1965) is characterized
  by low mortality and high fertility;  Madagascar (1966) is
  characterized by very high mortality and fertility, while the United
  States (1967) is characterized by low mortality and low fertility.
}
\source{
  Keyfitz and Flieger (1990), Goodman et al. (1974)
}
\references{
  Keyfitz, N., and W. Flieger. 1990. World population growth and
  aging: Demographic trends in the late twentieth century. Chicago:
  University of Chicago Press.

Goodman, L. A., N. Keyfitz, and T. W. Pullum. 1974. Family formation
and the frequency of various kinship relationships. Theoretical
Population Biology 5:1-27.
}

\examples{
data(goodman)
## Plot age-specific mortality rates for the three populations
with(goodman, plot(age, usa.nDx/usa.nKx, type="l", log="y", xlab="Age", ylab="nMx"))
with(goodman, lines(age, mad.nDx/mad.nKx, lty=2))
with(goodman, lines(age, ven.nDx/ven.nKx, lty=3))
legend("bottomright", c("USA","Madagascar","Venezuela"), lty=c(1,2,3))

}
\keyword{datasets}
