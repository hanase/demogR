`stoch.sens` <-
function(env,amat,k){

# A translation of Caswell's (2001) Matlab code fragment
  tlimit <- length(env)
  wvec <- rep(1/k,k)
  w <- cbind(wvec)

  # generate sequence of structure vectors

  r <- rep(0,tlimit)

  for(i in 1:tlimit){
    a <- amat[,env[i]]
    a <- matrix(a,nrow=k,ncol=k)
    wvec <- a%*%wvec
    r[i] <- sum(wvec)
    wvec <- wvec/r[i]
    w <- cbind(w,wvec)
  }

  # specifiy initial reproductive value vector

  vvec <- rep(1/k,k)
  v <- cbind(vvec)

  for(i in rev(1:tlimit)){
    a <- amat[,env[i]]
    a <- matrix(a,nrow=k,ncol=k)
    vvec <- vvec%*%a
    v <- cbind(t(vvec),v)
  }

  sensmat <- matrix(0,nrow=k,ncol=k)
  elasmat <- matrix(0,nrow=k,ncol=k)

  for(i in 1:tlimit){
    # for some reason, need the as.numeric() to get the division by
    # scalar to work 
    sensmat <- sensmat+((v[,i+1]%*%t(w[,i])) /
                        as.numeric(r[i]*t(v[,i+1])%*%w[,i+1]))
    a <- amat[,env[i]]
    a <- matrix(a,nrow=k,ncol=k)
    elasmat <- elasmat+((v[,i+1]%*%t(w[,i])*a) /
                        as.numeric((r[i]*t(v[,i+1])%*%w[,i+1])))
  }

  sensmat <- sensmat/tlimit
  elasmat <- elasmat/tlimit

  out <- list(sensitivities=sensmat, elasticities=elasmat)

  out
}
