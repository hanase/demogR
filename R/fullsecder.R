`fullsecder` <-
function(A){
  q <- A != 0
  size <- dim(A)
  qq <- matrix(q,nc=1)
  
  D <- NULL

  for(j in 1:size[2]){ # will work column-wise
    for(i in 1:size[1]){
      if(A[i,j]!=0){
        d2 <- secder(A,i,j)
        D <- cbind(D, matrix(d2,nc=1)*qq)
      }
    }
  }

  # find the non-zero elements
  qq <- which(D[,1] !=0)
  D <- D[qq,]
  dd <- dim(D)

  
  uu <- which(A>0, arr.ind=TRUE) # gives indices of non-zero elements
  o <- order(uu[,1])
  uu <- uu[o,]             # puts fertilities first
  D <- D[o,o]
  m <- length(uu[,1])
  # This makes dimnames for the final matrix to aid in interpretation 
  uuu <- rep(0,m)
  for(i in 1:m) uuu[i] <- paste(uu[i,1],uu[i,2],sep="")
  D <- matrix(D, nr=dd[1], nc=dd[2], dimnames=list(uuu,uuu))
  D
}

