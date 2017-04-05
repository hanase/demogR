\name{thar}
\alias{thar}
\docType{data}
\title{ thar }
\description{
  Caughley's (1966) survival data on Himalayan thar (\emph{Hemitargus jemlahicus}).
}
\usage{data(thar)}
\format{
  A data frame with 13 observations on the following 3 variables.
  \describe{
    \item{\code{age}}{ age class }
    \item{\code{count}}{ the number of individuals at the start of the
    age class}
    \item{\code{deaths}}{ observed deaths in the age class}
  }
}
\details{
 Caughley (1966) smoothed the counts within age classes and this data
 set used those reported smoothed counts.
}
\source{

Caughley, G. 1966. Mortality patterns in mammals. Ecology 47:906-918.


}

\examples{
data(thar)
thar.lt <- with(thar, 
    life.table(x=age, 
        nDx=deaths, 
        nKx=count, 
        type="cohort", 
        iwidth=1, 
        width12=c(1,1)))
thar.lt
}
\keyword{datasets}
