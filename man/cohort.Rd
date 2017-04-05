\name{cohort}
\alias{cohort}

\title{ cohort }
\description{
  Dummy function called in \code{life.table}.
}
\usage{
cohort(width12)
}

\arguments{
  \item{width12}{ width of the first two age classes }
}
\details{
  This is essentially a place-holder function used when passing the
  option \code{type="cohort"} to \code{life.table}.
}

\value{
A two element vector representing the first two values of the nax
column of the life table.
}

\seealso{  \code{\link{life.table}} }
\examples{
data(thar)
thar.lt <- with(thar, life.table(x=age, nDx=deaths, nKx=count, 
    type="cohort", iwidth=1, width12=c(1,1)))
thar.lt
}
\keyword{ array }

