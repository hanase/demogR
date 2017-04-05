\name{stoch.sens}
\alias{stoch.sens}

\title{ stoch.sens }
\description{
  Calculates the sensitivity of the stochastic growth rate to
  perturbations in the mean demographic projection matrix.
}
\usage{
stoch.sens(env, amat, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{ environmental sequence }
  \item{amat}{ matrices describing the population dynamics in the
  different environmental states, organized as columns of a matrix }
  \item{k}{ rank of the projection matrices }
}
\details{

See Caswell (2001, section 14.4) for details.

}
\value{
  A list with two elements:

  \item{sensitivities }{ sensitivities of the stochastic growth rate }
  \item{elasticities }{ elasticities of the stochastic growth rate }
  
}
\references{ 

Caswell, H. 2001. Matrix population models: Construction, analysis,
and interpretation. 2nd ed. Sunderland, MA: Sinauer. 

Haridas, C. V., and S. Tuljapurkar. 2005. Elasticities in variable
environments: Properties and implications. American Naturalist 166
(4):481-495. 

Tuljapurkar, S. 1990. Population dynamics in variable environments.
Vol. 85, Lecture notes in biomathematics. Berlin: Springer-Veralg.

Tuljapurkar, S., and C. V. Haridas. 2006. Temporal autocorrelation and
stochastic population growth. Ecology Letters 9 (3):324-334. 

}
\note{ 

This is essentially a translation of Caswell's (2001, section
14.4.1.1) MATLAB code.
}

\seealso{  \code{\link{lams}}, \code{\link{eigen.analysis}} }

\examples{

## Simulate an i.i.d. sequence of 3 environmental states

env <- floor(runif(100,0,3))+1
px1 <- rbeta(4,9.5,0.5)
px2 <- rbeta(4,7.5,2.5)
px3 <- rbeta(4,5,5)
mx <- c(0,rgamma(4,5,10))

A1 <- odiag(px1,-1)
A2 <- odiag(px2,-1)
A3 <- odiag(px3,-1)
A1[1,] <- leslie.row1(mx,px1)
A2[1,] <- leslie.row1(mx,px2)
A3[1,] <- leslie.row1(mx,px3)
amat <- cbind(matrix(A1,nr=25), matrix(A2,nr=25), matrix(A3,nr=25))
stoch.sens(env,amat,k=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
\keyword{ algebra }
\keyword{ ts }
