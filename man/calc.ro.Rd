\name{calc.ro}
\alias{calc.ro}

\title{ Caclulate net reproduction number from a demographic
projection matrix  }
\description{
  Calculate the net reproduction number ($R_0$) from an age or
  stage-classified demographic projection matrix.
}
\usage{
calc.ro(A, N.out = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ A demographic projection matrix }
  \item{N.out}{ Return the fundamental matrix (N) of the Markov chain? }
}
\details{
  Calculates the net reproduction number ($R_0$) from an age or
  stage-classified demographic projection matrix by first decomposing
  the k x k projection matrix A into two component matrices, T and F.  T
  collects the transitions between life-cycle stages while F collects
  the fertility transitions.  For an age-classified Leslie matrix, T
  will contain only the sub-diagonal of A and F will contain only the
  first row of A.  The fundamental matrix is given by $N = (I-T)^-1$,
  where I is a k x k identity matrix.  $R_0$ is the leading eigenvalue
  of the matrix FN.
}
\value{
If the (default) option N.out=FALSE is used, the net reproduction number is returned as a
single value.  If N.out=TRUE, the returned value is a list of two items:
  \item{ro }{Net reproduction number}
  \item{N }{Fundamental matrix}
  ...
}
\references{ 
Caswell, H. 2001. Matrix population models: construction, analysis,
and interpretation, Second edition. Sinauer, Sunderland,
Massachusetts, USA.
 }
\author{ James Holland Jones }

\seealso{ see \code{\link{leslie.matrix}} }
\examples{
## Create a Leslie matrix
Px <- c(0.77, 0.95, 0.98, 0.97)
Fx <- c(0,0,1,1.2,1)
L <- odiag(Px,-1)
L[1,] <- Fx
## Calculate $R_0$
calc.ro(L)
## look at the fundamental matrix
calc.ro(L, N.out=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }

