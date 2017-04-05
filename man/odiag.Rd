\name{odiag}
\alias{odiag}

\title{ odiag }
\description{
  Utility function that creates a square matrix with a vector on the off-diagonal or
  extracts the specified off-diagonal vector.
}
\usage{
odiag(A, at = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ either a vector or a matrix }
  \item{at}{ which diagonal? }
}
\details{
  If the first argument is a vector of length k, the result will be a
  square matrix of rank k+at with the provided vector along the 'at'
  diagonal.  Positive values for 'at' place the vector above the
  diagonal, negative values below the diagonal, and \code{at=0} (the
  default) places the vector on the diagonal.
}
\value{
A vector if argument \code{A} is a matrix and a matrix if \code{A} is
a vector.
}


\note{ 
\code{odiag} is meant to capture some of the functionality of the
MATLAB function diag().
}
\seealso{  \code{\link{diag}} }
\examples{

## Construct a matrix from a vector
## random survival probabilities with mean 0.9 and variance 0.0082

y <- rbeta(4,9,1)
A <- odiag(y,-1)

## add fertilities
F <- c(0,rep(1,4))
A[1,] <- F

## Extract a vector from a matrix
A <- matrix(rnorm(25), nr=5, nc=5)
odiag(A,2)
}
\keyword{ array }
\keyword{ algebra }
