\name{m2v}
\alias{m2v}

\title{ m2v }
\description{
  Utility that takes a matrix n x m and stacks its columns to form a
  column vector of length mn.
}
\usage{
m2v(A)
}

\arguments{
  \item{A}{ an m x n matrix}
}

\value{
 For an m x n input matrix, returns a m*n x 1 matrix (i.e., a column vector).

}

\note{ 

       Meant to capture the MATLAB functionality of the colon
       operator, A(:)
}
\seealso{  \code{\link{odiag}} }
\examples{
A <- matrix(rnorm(9), nr=3, nc=3)
m2v(A)
}
\keyword{ array }
\keyword{ algebra }
