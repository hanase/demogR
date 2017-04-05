`elassens` <-
function(A,k,l){
  rank <- dim(A)[1]
  ea <- eigen.analysis(A)
  lambda <- ea$lambda1
  s <- ea$sensitivities
  skl <- s[k,l]
  d2 <- secder(A,k,l)

  delta.kl <- matrix(0,nr=rank,nc=rank)
  delta.kl[k,l] <- 1

  es <- (A/lambda) * d2 - (A/lambda^2) * s*skl + delta.kl * s/lambda

  return(round(es,4))
}

