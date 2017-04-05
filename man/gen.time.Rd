\name{gen.time}
\alias{gen.time}

\title{ gen.time }
\description{
  Calculates the generation time for an age or stage-classified
  demographic projection matrix.
}
\usage{
gen.time(A, peryear = 5)
}

\arguments{
  \item{A}{ demographic projection matrix }
  \item{peryear}{ width of the age classes }
}
\details{
  Calculates the generation time (T) for an age or stage-classified
  demographic projection matrix using the identity 

  \eqn{R_0 = \exp(r*T)}

  where \eqn{R_0} is the net reproduction number and r is the intrinsic rate
  of increase (=log(lambda)).  

  Generation time is the amount of time that it takes a typical female
  to produce \eqn{R_0} offspring or, equivalently, the amount of time it
  takes a population growing with instantaneous rate r to increase by
  a factor of \eqn{R_0}.
}
\value{

The generation time implied by the demographic projection matrix.
}
\references{ 

Keyfitz, N., and H. Caswell. 2005. Applied mathematical
demography. 3rd ed. New York: Springer. 

Caswell, H. 2001. Matrix population models: Construction, analysis,
and interpretation. 2nd ed. Sunderland, MA: Sinauer. 

Preston, S.H., P. Heuveline, and F. Guillot. 2001. Demography:
 Measuring and modeling population processes. Oxford: Blackwell. 
 }

\note{ Calls function \code{calc.ro}, which calculates \eqn{R_0} from the
 fundamental matrix of the Markov transition matrix (Caswell 2001).

}
\seealso{  \code{\link{calc.ro}}, \code{\link{eigen.analysis}} }
\examples{

## compare generation times from two populations with very different
## R_0 values
data(goodman)
ven <- with(goodman, life.table(x=age, nKx=ven.nKx, nDx=ven.nDx))
ven.mx <- with(goodman, ven.bx/ven.nKx)
A <- leslie.matrix(lx=ven$nLx,mx=ven.mx)

usa <- with(goodman, life.table(x=age, nKx=usa.nKx, nDx=usa.nDx))
usa.mx <- with(goodman, usa.bx/usa.nKx)
B <- leslie.matrix(lx=usa$nLx,mx=usa.mx)

gen.time(A)
gen.time(B)
}
\keyword{ array }
\keyword{ algebra }
