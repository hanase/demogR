`leslie.row1` <-
function(mx,px, L=NULL, SRB=1.05, peryear=5, one.sex=TRUE){
 # takes the maternity function mx, and survivorship function
 # Lx and returns the fertilities for the first row of a Leslie 
 # matrix
 #
 # maternity function is yearly expected number of births

  if(is.null(L)) s <- sqrt(px[1]) * peryear
  else s <- L
  k <- length(px)+1
  Fx <- NULL
  for(i in 1:k-1){
     Fx[i] <- s * (mx[i] + px[i]*mx[i+1])/2
  }
  Fx[k] <- mx[k]/2
 if(one.sex) Fx <- Fx/(1+SRB)
Fx
}

