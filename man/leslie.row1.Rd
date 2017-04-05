\name{leslie.row1}
\alias{leslie.row1}

\title{ leslie.row1 }
\description{
 A utility for calculating the first row of a Leslie matrix from
 vectors of interval survival probabilities and age-specific fertility
 rates. 
}
\usage{
leslie.row1(mx, px, L = NULL, SRB = 1.05, peryear = 5, one.sex = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mx}{ vector of age-specific fertilities (length=k) }
  \item{px}{ vector of interval survival probabilities (length=k-1) }
  \item{L}{ person-years lived in the first interval; default
  \code{L=NULL} in which case each element is multiplied by
  peryear*sqrt(px[1]), otherwise each element is multiplied by L.  }
  \item{SRB}{ sex ratio at birth; default \code{SRB=1.05}. }
  \item{peryear}{ width of the age-class/projection interval }
  \item{one.sex}{ should each element be reduced by a factor
  1/(1+SRB)? default \code{one.sex=TRUE}. }
}
\details{
 Calculates the first-row entries for a Leslie matrix.  This utility
 is particularly useful for constructing simulated Leslie matrices. 
}
\value{
	A vector of length k, where k is the rank of the projection
	matrix A.
}

\seealso{  \code{\link{leslie.matrix}} }
\examples{
## simulate survival values from a beta(9,1) distribution

px <- rbeta(10,9,1)

## simulate fertility values from a gamma(2,10) distribution with age
##   at first reproduction = 10

mx <- c(0,0,0,rgamma(8,2,10))
A <- odiag(px,-1)
Fx <- leslie.row1(mx,px)
A[1,] <- Fx


}

\keyword{ array }
\keyword{ algebra }
