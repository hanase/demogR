`calc.ro` <-
function(A,N.out=FALSE){
  # Net reproduction number from Leslie matrix
  # assumes age-structured Leslie matrix
  k <- dim(A)[1]
  T <- A
  T[1,] <- 0                     # matrix of transitions
  F <- matrix(0,nrow=k,ncol=k)
  F[1,] <- A[1,]                 # matrix of births
  N <- solve(diag(k)-T)          # fundamental matrix
  ev <- eigen(F%*%N)
  imax <- which(ev$values==max(Re(ev$values)))
  ro <- ev$values[imax]          # same as FN[1,1]

  if(N.out) out <- list(ro, N)
  else out <- ro

  out
}

