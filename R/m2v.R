`m2v` <-
function(A){
  k <- dim(A)
  k1 <- k[1]
  k2 <- k[2]
  s <- matrix(A,nr=k1*k2,nc=1)
  s
}

