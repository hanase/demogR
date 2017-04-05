`lams` <-
function(aseq,n=5) {
  k <- sqrt(dim(aseq)[1])
  abar <- apply(aseq,1,mean)
  abar <- matrix(abar,nrow=k,ncol=k)
  ea <- eigen.analysis(abar)
  lambda <- ea$lambda1
  s <- ea$sensitivities
  s <- matrix(s,nrow=k*k)

  ## covariace matrix

  c <- cov(t(aseq))
  lams <- log(lambda) - (t(s) %*% c %*% s)/(2*lambda^2)
  lams <- lams/n

  lams
}

