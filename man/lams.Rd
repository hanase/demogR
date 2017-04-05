\name{lams}
\alias{lams}

\title{ lams }
\description{
  Calculates the logarithm of the stochastic growth rate using
  Tuljapurkar's second order approximation for independent and
  identically distributed environments.
}
\usage{
lams(aseq, n = 5)
}

\arguments{
  \item{aseq}{ sequence of matrices with each matrix given as a
  re-shaped column of aseq }
  \item{n}{ width of the projection interval/age-class }
}
\details{
  Uses Tuljapurkar's second order approximation for independent and
  identically distributed (i.i.d.) environments.
}
\value{
  The long-run growth rate for the population with projection matrices
  given by aseq.
}
\references{ 
Tuljapurkar, S. 1990. Population dynamics in variable
  environments. Edited by S. A. Levin. Vol. 85, Lecture notes in
  biomathematics. Berlin: Springer-Veralg. 

Caswell, H. 2001. Matrix population models: Construction, analysis,
and interpretation. 2nd ed. Sunderland, MA: Sinauer. 
}

\seealso{  \code{\link{stoch.sens}}, \code{\link{eigen.analysis}} }
\examples{
## simulate two matrices: good year/bad year
## beta distributed survival, gamma fertility

px1 <- rbeta(5,9,1)
px2 <- rbeta(5,7,3)
mx1 <- c(0,rgamma(5,10,10))
mx2 <- c(0,rgamma(5,7,10))
## good year matrix
A1 <- odiag(px1,-1)
A1[1,] <- mx1

## bad year matrix
A2 <- odiag(px2,-1)
A2[1,] <- mx2
aseq <- cbind(matrix(A1,nr=36,nc=1), matrix(A2,nr=36,nc=1))
lams(aseq)
}


\keyword{ array}
\keyword{ algebra }
