\name{coale}
\alias{coale}

\title{ coale }
\description{
Utility used by \code{life.table} to graduate the person-years lived
  by those dying in the interval by the method of Coale and Demeny. 
}
\usage{
coale(b1, b4, nMx)
}

\arguments{
  \item{b1}{ two element vector of regression coefficients for
  graduating 1a0 provided in  Coale et al. (1983) }
  \item{b4}{ two element vector of regression coefficients for
  graduating 4a1 provided in  Coale et al. (1983) }
  \item{nMx}{ period central death rates: \code{=nDx/nKx} }
}
\details{
  Utility function called by \code{life.table}.
}
\value{
The first two values (age classes 0-1 and 1-5) of the nax column of
  a period life table.
}
\references{ 

Coale, A., P. Demeny, and B. Vaughn. 1983. Regional model life tables
and stable populations. 2nd ed. New York: Academic Press.

 }

\seealso{  \code{\link{life.table}} }
\examples{
data(goodman)
vlt1 <- with(goodman, life.table(x=age, nKx=ven.nKx, nDx=ven.nDx, type="cd"))

}

\keyword{ array}

