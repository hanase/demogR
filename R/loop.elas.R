`loop.elas` <-
function(A,draw.plot=TRUE,peryear=5, xlab="Loop Elasticity",
                      ylab="Age", xlim=c(0,(maxe+0.02)), ...){
  k <- dim(A)[1]
  ea <- eigen.analysis(A)
  e <- ea$elasticities
  ef <- e[1,]

  loop.e <- ef*1:k

  maxe <- max(loop.e)

  if(draw.plot){
    barplot(loop.e,
            names.arg=peryear*1:k,
            horiz=TRUE,
            beside=TRUE,
            xlim=xlim,
            xlab=xlab,
            ylab=ylab,
            ...)
  }

  return(loop.e)
}

