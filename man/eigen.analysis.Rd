\name{eigen.analysis}
\alias{eigen.analysis}

\title{ eigen.analysis }
\description{
  Calculates the asymptotic growth rate and related quantities from a
  demographic projection matrix.
}
\usage{
eigen.analysis(A)
}

\arguments{
  \item{A}{ demographic projection matrix }
}
\details{

 Calculates the asymptotic growth rate (lambda) of a population
  described by demographic projection matrix A.  The asymptotic growth
  rate of the population is given by the dominant eigenvalue of the
  projection matrix.  By the Perron-Frobenius Theorem, this eigenvalue
  is guaranteed to be real, positive and strictly greater than all the
  other eigenvalues if the matrix A is non-negative, irreducible, and
  primitive (for details see Caswell (2001)).

  Also calculates the damping ratio (rho), eigenvalue sensitivities,
  eigenvalue elasticities, the stable age distribution (for the
  communicating parts of the life cycle), and  scaled reproductive
  values. 

  The damping ratio is the ratio of the dominant eigenvalue and the
  absolute value of the second eigenvalue.  rho is a measure of the
  rate of convergence to the stable age-distribution.  A population
  characterized by damping ratio rho will converge asymptotically to
  the stable age distribution exponentially with rate at least as fast
  as log(rho).  Clearly, a population already at or very near the
  stable age distribution will converge faster, but rho provides an
  upper bound.

  The eigenvalue sensitivities are the partial derivatives of lambda
  with respect to a perturbation in matrix element $a_ij$.  The
  sensitivities measure the selection gradient on the life-cycle
  (Lande 1982).   The eigenvalue elasticities are scaled to be
  proportional sensitivities of lambda to a perturbation in $a_ij$.
  Elasticities have a number of desirable properties including, their
  sum across all life-cycle transitions is unity and the sum of the
  elasticities of all  incoming arcs to a life-cycle stage must equal
  the sum of all outgoing arcs (van Groenendael et al 1994).

   The stable age distribution is normalized to represent the
   proportion in each of the communicating age classes.  If the
   population is characterized by post-reproductive survival (and
   hence age classes that do not communicate with the rest of the life
   cycle graph), then other methods should be used to calculate to
   stable distribution.  For example, from classic stable population
   theory, we know that the stable age distribution of the population
   c(x) is given by the relationship: 

   c(x) = b l(x) exp(-r*x)

   where b is the gross birth rate, l(x) is survivorship to age x and
   r is the rate of increase of the population (=log(lambda)).  See
   Coale (1972) or Preston et al. (2001) for details.

   The age-specific reproductive values are normalized so that the
   reproductive value of the first age class is unity.  Problems
   associated with post-reproductive survival are irrelevant for
   reproductive value since the reproductive value of
   post-reproductive individuals is, by definition, zero. 
}
\value{
  A list with six components:

  \item{lambda1 }{ the asymptotic growth rate (dominant eigenvalue) of
  A }
  \item{rho }{ damping ratio of A }
  \item{sensitivities}{ eigenvalue sensitivities of A }
  \item{elasticities}{ eigenvalue elasticities of A }
  \item{stable.age}{  stable age distribution of A }
  \item{repro.value}{  reproductive values of A }

}
\references{ 

Caswell, H. 2001. Matrix population models: Construction, analysis,
and interpretation. 2nd ed. Sunderland, MA: Sinauer. 

Coale, A.J. 1972. The growth and structure of human populations: A
mathematical investigation. Princeton: Princeton University Press. 

Lande, R. A. 1982. A quantitative genetic theory of life history
evolution. Ecology 63:607-615. 

van Groenendael, J., H. De Kroon, S. Kalisz, and
S. Tuljapurkar. 1994. Loop analysis: Evaluating life history pathways
in population projection matrices. Ecology 75 (8):2410-2415. 
}

\seealso{  \code{\link{leslie.matrix}}, \code{\link{secder}} }
\examples{

data(goodman)
ven <- with(goodman, life.table(x=age, nKx=ven.nKx, nDx=ven.nDx))
ven.mx <- with(goodman, ven.bx/ven.nKx)
A <- leslie.matrix(lx=ven$nLx,mx=ven.mx)

usa <- with(goodman, life.table(x=age, nKx=usa.nKx, nDx=usa.nDx))
usa.mx <- with(goodman, usa.bx/usa.nKx)
B <- leslie.matrix(lx=usa$nLx,mx=usa.mx)

eav <- eigen.analysis(A)
eau <- eigen.analysis(B)

## compare rates of increase

eav$lambda1
eau$lambda1

## compare stable age distributions for Venezuela and the USA
op <- par(no.readonly = TRUE)
age <- seq(0,50,by=5)
par(mfrow=c(1,2))
barplot(eav$stable.age, 
             names.arg=age, 
	     horiz=TRUE, 
	     col="grey", 
	     ylab="Age")
barplot(eau$stable.age, 
	     names.arg=age[-11], 
	     horiz=TRUE, 
	     col="blue", 
	     ylab="Age", 
	     xlim=c(0,0.2))

par(op)

## Compare reproductive values for Venezuela and the USA

plot(age, eav$repro.value, type="l", xlab="Age", ylab="Reproductive Value")
lines(age[-11], eau$repro.value, lty=2)
legend("topright", c("Venezuela (1965)", "USA (1967)"), lty=c(1,2))
}

\keyword{ array  }
\keyword{ algebra }
