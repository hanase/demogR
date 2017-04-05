`keyfitz` <-
function(b0,nMx){
  nax12 <- c(0,0)
  nax12[1] <- b0[1] + b0[2] *nMx[1]
  nax12[2] <- 1.5
  return(nax12)
}

