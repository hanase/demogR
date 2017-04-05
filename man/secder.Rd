\name{secder}
\alias{secder}

\title{secder }
\description{
  Calculates the second derivatives of the dominant eigenvalue of the
  demographic projection matrix for all non-zero transitions with
  respect to one specified transition.
}
\usage{
secder(A, k, l)
}

\arguments{
  \item{A}{ demographic projection matrix }
  \item{k}{ row index for the specified transition }
  \item{l}{ column index for the specified transition }
}
\details{
   See Caswell (1996, 2001) for details on second derivatives of the
  dominant eigenvalue.
}
\value{

  A square matrix of the same rank as A where each element $s_ij$ is the
  second derivative of the dominant eigenvalue of A, $\partial^2
  \lambda/\partial a_{ij} \partial a_{kl}$.

}
\references{  

Caswell, H. 1996. Second derivatives of population growth rate:
Calculation and applications. Ecology 77 (3):870-879.

Caswell, H. 2001. Matrix population models: Construction, analysis,
and interpretation. 2nd ed. Sunderland, MA: Sinauer.

 }

\note{

The eigenvalue second derivatives are essential for calculating both
perturbation analyses of the eigenvalue elasticities and stochastic
sensitivities.  \code{secder} is used in functions to calculate both
these quantities.

}

\seealso{  \code{\link{fullsecder}}, \code{\link{elassens}},
\code{\link{eigen.analysis}}, \code{\link{stoch.sens}} } 
\examples{

## eigenvalue second derivatives of the US projection matrix from 1967
## with respect to infant survival
data(goodman)
ult <- with(goodman, life.table(x=age, nKx=usa.nKx, nDx=usa.nDx))
mx <- goodman$usa.bx/goodman$usa.nKx
usa <- leslie.matrix(lx=ult$nLx,mx=mx)

sd21 <- secder(usa,2,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
\keyword{ algebra }
