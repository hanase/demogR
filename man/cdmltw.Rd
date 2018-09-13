\name{cdmltw}
\alias{cdmltw}
\alias{cdmlte}
\alias{cdmlts}
\alias{cdmltn}

\title{ Coale-Demeny Model Life Tables }
\description{
  Coale-Demeny regional model life tables
}
\usage{
cdmltw(sex = "F")
cdmlte(sex = "F")
cdmlts(sex = "F")
cdmltn(sex = "F")
}

\arguments{
  \item{sex}{ female or male life tables.  Defaults to 'F' (female). }
}

\details{

  Provides region-specific model life tables for 25 different levels of mortality,
  indexed by life expectancy at age 10 (with 1 being the lowest and 25
  being the highest life expectancy at age 10). 
  
  The regions are West (\code{cdmltw}), East (\code{cdmlte}), South (\code{cdmlts}), and North (\code{cdmltn}).

}

\value{

  A list containing two 21 element vectors, one 25 element vector and eight 25 x 21 matrices:
  25 mortality levels by 21 age classes.  The eight matrices contain
  the standard columns of a life table for each mortality level:

  \item{ age }{ the 21 age classes }
  \item{ width }{widths of the 21 age classes}
  \item{ e10 }{ the 25 mortality levels as life expectancy at age 10}
  \item{ lx }{survival probabilities to exact age x }
  \item{ nqx }{ probabilities of death between ages x and x+n}
  \item{ nax }{ person-years lived by those dying between ages x and
  x+n }
  \item{ ndx }{ proportion of all deaths occurring between ages x and
  x+n }
  \item{ nLx }{ person-years lived in the interval x to x+n }
  \item{ nmx }{ central death rate in the interval x to x+n}
  \item{ Tx }{ person-years of life remaining to the cohort at age x }
  \item{ ex }{ life expectancy at age x }
}
\references{ 

Coale, A., P. Demeny, and B. Vaughn. 1983. Regional model life tables
and stable populations. 2nd ed. New York: Academic Press.

 }
\author{ Kenneth Wachter, Department of Demography, University of
California, Berkeley.

Modified by James Holland Jones, Department of Anthropological
Sciences, Stanford University

Corrected by Jim Oeppen, Southern Denmark University, Odense, Denmark
 }

\seealso{  \code{\link{life.table}} }
\examples{
## Generate West model life table family for females
Wf <- cdmltw()

## Plot survivorship for the 25 levels
plot(Wf$age,Wf$lx[1,], type="l", xlab="Age", ylab="Survivorship")
for(i in 2:25) lines(Wf$age, Wf$lx[i,])
}
\keyword{ array }
\keyword{ algebra }
