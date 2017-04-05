\name{project.leslie}
\alias{project.leslie}

\title{ project.leslie }
\description{
  Projects an population vector tmax intervals by
  pre-multiplication with a Leslie matrix.
}
\usage{
project.leslie(A, no, tmax, pop.sum = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ a k x k projection matrix }
  \item{no}{ a k x 1 population vector }
  \item{tmax}{ number of time steps to project the vector }
  \item{pop.sum}{ logical; If 'TRUE', the age-classes of the projected
  population are summed, yielding a single total population vector}
}
\details{
  Takes an initial population vector, no, and pre-multiplies by the
  demographic projection matrix, A, tmax times.  This projection will
  be tmax*n years into the future, where n is the width of the
  age-classes in the Leslie matrix, A.
}

\value{

	If \code{pop.sum=FALSE} (the default), the value will be a k x
	tmax+1 matrix. The first column of the matrix is no and each
	subsequent column represents the population structure at time
	step 1, 2, ..., tmax.

	If  \code{pop.sum=TRUE}, the value will be a vector of length
	tmax+1, where each element of the vector is the total
	population at time t=0, 1, ..., tmax.

}

\seealso{  \code{\link{leslie.matrix}}, \code{\link{eigen.analysis}} }
\examples{
data(goodman)
ult <- with(goodman, life.table(x=age, nKx=usa.nKx, nDx=usa.nDx))
mx <- goodman$usa.bx/goodman$usa.nKx
usa <- leslie.matrix(lx=ult$nLx,mx=mx)

## initial population is the age structure in 1967
## this includes the number of 0-1 year-olds and the number of 1-4
## year-olds, which need to be combined into a single 5-year wide
## element 

no <- goodman$usa.nKx[3:11]
no <- c(sum(goodman$usa.nKx[1:2]),no)/1e6

N <- project.leslie(usa,no,tmax=20)
plot( seq(1967,2067,by=5), apply(N,2,sum), 
         type="l", 
	 xlab="Years", 
	 ylab="Projected Population (millions)")
title("Projected Total Population Size (Ages 0-45) of the United States")

##   Note that this is a lousy projection since it assumes constant
##   demographic rates and a closed population and only accounts for
##   ages 0-45
}

\keyword{ array }
\keyword{ algebra }
