\name{life.table}
\alias{life.table}

\title{ life.table }
\description{
  Construct either a period or cohort life table from enumerated deaths
  and mid-interval population estimates.
}
\usage{
life.table(x, nDx, nKx, 
	      b0 = c(0.07, 1.7), 
	      b1 = c(0.053, 2.8), 
	      b4 = c(1.522, 1.518), 
	      type = "kf", 
	      nxx = 0, 
	      iwidth = 5, 
	      width12 = c(1, 4))
}

\arguments{
  \item{x}{ age at the beginning of the age classes of the life table }
  \item{nDx}{ deaths }
  \item{nKx}{ population size }
  \item{b0}{ coefficients used in Keyfitz-Flieger graduation }
  \item{b1}{ first set of coefficients used in Coale-Demeny graduation }
  \item{b4}{ second set of coefficients used in Coale-Demeny graduation }
  \item{type}{ type of life table calculation: "kf", "cd", or
  "cohort". Default is "kf". }
  \item{nxx}{ person-years lived by those dying in the last (possibly
  open) age-class. If \code{nxx=0}, the person-years lived by those
  dying in the interval is the inverse of the central death rate
  (corresponding to exponentially distributed failure times). }
  \item{iwidth}{ width of the age intervals. Default \code{iwidth=5}. }
  \item{width12}{ width of the first two age classes. Default
  \code{width12=c(1,4)}.  }
}
\details{

Constructs a period or cohort life tables from enumerated deaths and
mid-interval population sizes (period) or enumerated deaths and
person-years at risk (cohort).  \code{x}, \code{nDx}, and \code{nKx} must all
the be same length.

There are currently three options for life table construction.  The
first two are for the construction of period life tables.  They differ
only in the way that person-years lived by those dying in the first
two intervals are handled.  For \code{type="kf"}, the default, the
first two values of nax estimated using Keyfitz and Fleiger's (1990)
regression method. For \code{type="cd"}, Coale and Demeny's method
(1983) is used.  The Coale-Demeny method uses different coefficients
depending on the level of early mortality.  As a result, this method
may work better for high-mortality populations. 

The third type of life table is a cohort life table, for which the
conversion from mortality rates to probabilities of death is
unnecessary, so the nax column of the life table is of limited
interest. 

}
\value{
 A dataframe with nine columns:

  \item{ x }{ age at the beginning of the interval }
  \item{ nax }{ person-years lived by those dying in the interval x to
  x+n }
  \item{nMx}{ period central death rate }
  \item{ nqx}{ probability of death between ages x and x+n }
  \item{ lx }{ probability of survival to exact age x }
  \item{ ndx }{ proportion of deaths occurring between ages x and x+n
  } 
  \item{ nLx }{ person-years lived in the interval x to x+n }
  \item{ Tx }{ person-years of life left in the cohort at age x }
  \item{ ex }{ life expectancy at age x}

}
\references{ 
Keyfitz, N. 1977. Introduction to the mathematics of populations. 2nd
ed. Menlo Park: Addison-Wesley.

Coale, A., P. Demeny, and B. Vaughn. 1983. Regional model life tables
and stable populations. 2nd ed. New York: Academic Press.

Keyfitz, N., and W. Flieger. 1990. World population growth and
  aging: Demographic trends in the late twentieth century. Chicago:
  University of Chicago Press.

Preston, S.H., P. Heuveline, and F. Guillot. 2001. Demography:
Measuring and modeling population processes. Oxford: Blackwell. 
 }


\seealso{  \code{\link{cdmltw}} }

\examples{
data(goodman)
## default type="kf"
vlt <- with(goodman, life.table(x=age, nKx=ven.nKx, nDx=ven.nDx))

## compare nax values for cd vs kf life tables
vlt1 <- with(goodman, life.table(x=age, nKx=ven.nKx, nDx=ven.nDx, type="cd"))
vlt$nax[1:2]
vlt1$nax[1:2]

## now construct a cohort life table for Himalayan thar,
##  (Hemitargus jemlahicus) 
data(thar)
thar.lt <- with(thar, life.table(x=age, nDx=deaths, nKx=count, 
    type="cohort", iwidth=1, width12=c(1,1)))
thar.lt
}

\keyword{ array }

