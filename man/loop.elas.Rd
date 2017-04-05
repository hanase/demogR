\name{loop.elas}
\alias{loop.elas}

\title{ loop.elas }
\description{
  Calculate  and plot loop elasticities for a demographic projection matrix .
}
\usage{
loop.elas(A, 
	       draw.plot = TRUE, 
	       peryear = 5, 
	       xlab = "Loop Elasticity", 
	       ylab = "Age", 
	       xlim = c(0, (maxe + 0.02)), 
	       ...)
}

\arguments{
  \item{A}{ demographic projection matrix }
  \item{draw.plot}{ should a plot be drawn? default:
  \code{draw.plot=TRUE} } 
  \item{peryear}{ width of the age-class/projection interval }
  \item{xlab}{ x-axis label; default \code{ylab="Age"}}
  \item{ylab}{ y-axis label; default \code{ylab="Loop Elasticity"}}
  \item{xlim}{ plotting limits for x-axis }
  \item{\dots}{ arguments to be passed to/from other methods }
}
\details{
  
 Calculates and plots the life-cycle decomposition of van Groenendael
 et al (1994). 

}

\value{
  A vector of length k (the rank of the projection matrix A) giving
  the loop elasticities for each of the loops in the life-cycle
  graph. 

}
\references{ 

Van Groenendael, J., H. De Kroon, S. Kalisz, and
S. Tuljapurkar. 1994. Loop analysis: Evaluating life history pathways
in population projection matrices. Ecology 75 (8):2410-2415. 


}

\seealso{  \code{\link{eigen.analysis}}, \code{\link{leslie.matrix}} }
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

## Compare the loop elasticities of Venezuela (1965) and the USA (1967)
le.usa <- loop.elas(A)
le.usa <- loop.elas(B)


}

\keyword{ array }
\keyword{ algebra }
