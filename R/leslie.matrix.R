`leslie.matrix` <-
function(lx,mx, L=TRUE, peryear=5, one.sex=TRUE,
                          SRB=1.05, infant.class=TRUE) {
# assumes lx and mx same length

# lx for quinquennia (as typical for humans) will typically be nLx --
# person-years lived in the 5-year interval; mx is the births/exposure;
# peryear=1 means that age-specific fertilities have been
# appropriately cumulated for the interval, if peryear=5, mx values
# will be multiplied by 5, etc.  

  # L=FALSE indicates that lx is used.  If L=TRUE, then nLx is used.
  # In the first case, lx[1] needs to be removed before calculating Px
  # while in the second, nLx[1] and nLx[2] need to be combined

  len1 <- length(lx)
  len2 <- length(mx)
  if(len1>len2) {
    warning("length of lx greater than the length of mx,\n lx truncated to length of mx")
    lx <- lx[1:len2]
  }
  if(len2>len1) {
    mx <- mx[1:len1]
  }


  
  if(infant.class)   mx <- mx[-2]
  fages <- which(mx>0)
  k <- max(fages)
  mx <- mx[1:k]


  # make the first-row entries for the Leslie matrix
  if(L){
    L1 <- lx[1]+lx[2]
    s <- L1
    lx <- c(L1,lx[-c(1,2)])
  }
   else{
     lx <- lx[-2]
     s <- sqrt(lx[2])*peryear
   }
                                       
  # calculate px

  px <- exp(diff(log(lx)))
  px <- px[1:(k-1)]

  # calculate Fx 
  
  Fx <- NULL
  for(i in 1:k-1) {
      Fx[i] <- s * (mx[i] + px[i]*mx[i+1])/2
  }
  Fx <- c(Fx,mx[k])

  # one-sex population
  if(one.sex)  Fx <- Fx/(1+SRB)
  A <- matrix(0, nrow=k, ncol=k)
  A[row(A) == col(A)+1] <- px
  A[1,] <- Fx

  class(A) <- "leslie.matrix"
  
  A
}

