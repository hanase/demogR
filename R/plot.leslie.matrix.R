`plot.leslie.matrix` <-
  function(x, y=NULL,
                        ..., 
                        main="", sub="",
                        ann=TRUE,
                        xlim=NULL, ylim=NULL,
                        axes=TRUE, 
                        col=c("black","grey"),
                        lwd=2,
                        xlab="Age", ylab="Sensitivity",
                        peryear=5) {
  A <- x                                ## just for consistency with other usage
  d <- dim(A)[1]
  ss <- c(A[row(A) == col(A)+1], NA)    ## add NA 
  fs <- A[1,]
  age <- peryear*(0:(d-1))
  isz <- fs==0        # determine non-zero fertilities
  
  xy <- xy.coords(age,ss,xlab=xlab,ylab=ylab)
  if(is.null(xlim))
    xlim <- range(xy$x[is.finite(xy$x)])
  if(is.null(ylim))
    ylim <- range(xy$y[is.finite(xy$y)])
  
  plot.new()
  plot.window(xlim,ylim, ...)
  lines(xy$x, xy$y, lwd=lwd, col=col[1], ...)
  lines(age[!isz], fs[!isz], lwd=lwd, col=col[2])
  if(axes) {
    axis(1)
    axis(2)
    box()
  }
  if (ann)
    title(main=main, sub=sub, xlab=xy$xlab, ylab=xy$ylab, ...)
}


