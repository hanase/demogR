`life.table` <-
function(x, nDx, nKx, b0=c(0.07,1.7), b1=c(0.053,2.8),
  b4=c(1.522,1.518), type="kf",nxx=0, iwidth=5, width12=c(1,4)){

 nmax <- length(x)

# period rates n_M_x

 nMx <- nDx/nKx

# n_a_x


 nax <- NULL

 nax12 <- switch(type,
                cd = coale(b1,b4,nMx),
                kf = keyfitz(b0,nMx),
                cohort = cohort(width12)
                 )


 nax[1] <- nax12[1]
 nax[2] <- nax12[2]
 

# nxx closes out the life table.  For low-mortality populations (e.g.,
# where there is under-reporting of old-age deaths), it is essential
# that nxx be specified.
   
 nax[3:(nmax-1)] <- iwidth/2
 if(nxx == 0)  nax[nmax] <- 1/nMx[nmax]
 else nax[nmax] <- nxx

 n <- c(width12, rep(iwidth, nmax - 3),999) #width of the intervals
 nqx <- (n*nMx) / (1 + (n-nax)*nMx)
 nqx[nmax] <- 1.0


 # preceding is unnecessary for type="cohort"
 if(type=="cohort")   nqx <- nMx
 
# survivorship lx

 lx <- cumprod(c(1,1-nqx))

 # l_{x+n}
 lxpn <- lx[-1]
 ndx <- -diff(lx)

# person-years lived by survivors
 nLx <- n*lxpn + ndx*nax
 Tx <- rev(cumsum(rev(nLx)))

 ex <- Tx/lx[1:nmax]

# format for printing

 lt <- data.frame(x, nax = round(nax,4), 
                nMx = round(nMx,4), 
                nqx = round(nqx[1:nmax],4), 
                lx = round(lx[1:nmax],4), 
                ndx = round(ndx,4), 
                nLx = round(nLx,4),
                Tx = round(Tx,2),
                ex = round(ex,2) )

lt
}

