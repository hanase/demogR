\name{elassens}
\alias{elassens}

\title{elassens }
\description{
  Calculates the partial derivatives of the eigenvalue elasticities of
  a demographic projection matrix.
}
\usage{
elassens(A, k, l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ a demographic projection matrix }
  \item{k}{ row index of the transition that the partial is calculated
  with respect to}
  \item{l}{ column index of the transition that the partial is calculated
  with respect to}
}
\details{
  See Caswell (2001) section 9.7.1.
}
\value{
A matrix of class "leslie.matrix".
}
\references{ 

Caswell, H. 2001. Matrix population models: Construction, analysis, and interpretation. 2nd ed. Sunderland, MA: Sinauer.

 }

\seealso{  \code{\link{secder}}, \code{\link{leslie.matrix}},
\code{\link{eigen.analysis}} } 
\examples{
data(goodman)
mlt <- with(goodman, life.table(x=age, nDx=mad.nDx, nKx=mad.nKx))
mx <- goodman$mad.bx/goodman$mad.nKx

##  make the Leslie matrix
mad <- leslie.matrix(lx=mlt$lx, mx=mx)

## now calculate the sensitivities of the elasticitities with respect
##to infant survival

es21 <- elassens(mad,2,1)
}
\keyword{array }
\keyword{ algebra }
