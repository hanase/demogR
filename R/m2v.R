`m2v` <-
function(A){
  k <- dim(A)
  k1 <- k[1]
  k2 <- k[2]
  s <- matrix(A,nrow=k1*k2,ncol=1)
  s
}

