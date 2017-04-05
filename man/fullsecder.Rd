\name{fullsecder}
\alias{fullsecder}

\title{ fullsecder }
\description{
  Calculates the second derivatives of the dominant eigenvalue of a
  square matrix, A, with respect to all non-zero elements of A.
}
\usage{
fullsecder(A)
}

\arguments{
  \item{A}{ a demographic projection matrix }
}
\details{
  Produces a matrix the columns of which are re-shaped matrices of
  second derivatives of the dominant eigenvalue of the projection
  matrix A with respect to each non-zero element in A.

  See Caswell (1996, 2001) for details on second derivatives of the
  dominant eigenvalue.

}
\value{
 A square matrix.  If A is a Leslie matrix of rank k, then the maximum
  rank of the resulting matrix is 2k-1 (since there are at most k
  non-zero fertilities and k-1 survival probabilities in a k x k
  Leslie matrix). 
}
\references{ 

Caswell, H. 1996. Second derivatives of population growth rate:
Calculation and applications. Ecology 77 (3):870-879.

Caswell, H. 2001. Matrix population models: Construction, analysis,
and interpretation. 2nd ed. Sunderland, MA: Sinauer.

}


\seealso{  \code{\link{secder}}, \code{\link{elassens}},
\code{\link{eigen.analysis}}, \code{\link{stoch.sens}}  } 
\examples{

## eigenvalue second derivatives of the US projection matrix from 1967
data(goodman)
ult <- with(goodman, life.table(x=age, nKx=usa.nKx, nDx=usa.nDx))
mx <- goodman$usa.bx/goodman$usa.nKx
usa <- leslie.matrix(lx=ult$nLx,mx=mx)

fs <- fullsecder(usa)

## plot the survival cross-second derivatives of lambda with respect
## to infant survival

plot( seq(0,40,by=5), fs["21",10:18], type="l", 
	 xlab="Age (j)", 
	 ylab=expression(paste(partialdiff^2 , lambda, "/", 
	 partialdiff, P[1] , partialdiff, P[j])))
abline(h=0,lty=3)
}
\keyword{ array }
\keyword{ algebra }
