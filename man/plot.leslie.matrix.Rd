\name{plot.leslie.matrix}
\alias{plot.leslie.matrix}
\title{ plot.leslie.matrix }
\description{
  Plots Leslie matrices and other objects of class "leslie.matrix".
}
\usage{
\method{plot}{leslie.matrix}(x, y=NULL,
                  ...,
                  main="", sub="",
                  ann=TRUE,
                  xlim=NULL, ylim=NULL,
                  axes=TRUE,
                  col=c("black","grey"),
                  lwd=2,
                  xlab="Age", ylab="Sensitivity",
                  peryear=5 )
}

\arguments{
  \item{x}{ demographic projection matrix }
  \item{y}{ y axis values; NULL}
  \item{\dots}{ arguments to be passed to/from other methods }
  \item{main}{ an overall title for the plot }
  \item{sub}{ a subtitle for the plot }
  \item{ann}{ logical; if TRUE, annotate plot  }
  \item{xlim}{plotting limits for x-axis; if \code{xlim=NULL}, range of
    x values used }
  \item{ylim}{ plotting limits for y-axis; if \code{xlim=NULL}, range of
    y values used }
  \item{axes}{logical; if TRUE, draw axes and box}
  \item{col}{ vector (length 2) of line colors for the plot; default
    \code{c("black","grey")}}
  \item{lwd}{ thickness of the plotted lines; default \code{lwd=2}}
  \item{peryear}{ width of the age-class/projection interval; default \code{peryear=5} }
  \item{xlab}{ x-axis label; default \code{xlab="Age"}}
  \item{ylab}{ y-axis label; default \code{xlab="Sensitivity"}}
}
\details{
  Plots the subdiagonal and the first row of an object of class
  leslie.matrix on common axes.  Probably most useful for plotting
  either eigenvalue sensitivities or elasticities of a projection matrix.
}
\value{
  NULL

}

\seealso{  \code{\link{eigen.analysis}}, \code{\link{loop.elas}} }
\examples{
data(goodman)
ven <- with(goodman, life.table(x=age, nKx=ven.nKx, nDx=ven.nDx))
ven.mx <- with(goodman, ven.bx/ven.nKx)
A <- leslie.matrix(lx=ven$nLx,mx=ven.mx)
ea.ven <- eigen.analysis(A)

usa <- with(goodman, life.table(x=age, nKx=usa.nKx, nDx=usa.nDx))
usa.mx <- with(goodman, usa.bx/usa.nKx)
B <- leslie.matrix(lx=usa$nLx,mx=usa.mx)
ea.usa <- eigen.analysis(B)

## Compare sensitivities for Venezuela and the USA

plot(ea.ven$sens, main="Venezuela")
plot(ea.usa$sens, main="USA")
}
\keyword{ array }
\keyword{ algebra }
