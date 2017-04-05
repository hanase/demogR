\name{keyfitz}
\alias{keyfitz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ keyfitz }
\description{
  Utility used by \code{life.table} to graduate the person-years lived
  by those dying in the interval by the method of Keyfitz and Flieger
  (1990). 
}
\usage{
keyfitz(b0, nMx)
}

\arguments{
  \item{b0}{ two element vector of regression coefficients provided in
  Keyfitz and Flieger (1990).  Default value: \code{b0=c(0.07,1.7)}. }
  \item{nMx}{ period central death rates: \code{=nDx/nKx}. }
}
\details{
  Utility function called by \code{life.table}.
}
\value{
 The first two values (age classes 0-1 and 1-5) of the nax column of
  a period life table.

}

\references{ 
Keyfitz, N., and W. Flieger. 1990. World population growth and
  aging: Demographic trends in the late twentieth century. Chicago:
  University of Chicago Press.
 }



}
\seealso{  \code{\link{life.table}} }
\examples{
data(goodman)
## default type="kf"
vlt <- with(goodman, life.table(x=age, nKx=ven.nKx, nDx=ven.nDx))


}
\keyword{ array }

