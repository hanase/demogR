`cdmltn` <-
function(sex="F"){
  if (sex != "F" & sex !="M"){
    stop("sex must be either F or M!")
  }
  
  # CD lifetables indexed by e_10 = eten
  # must do a first pass of females to get eten values for both sexes
  for(i in 1:2){
    if(i==1){
       eten <- seq(20,84, by = 2)
     }else eten <- etenf



#    Level 1 is ef=20,  Level 25 is  ef = 80
#    Estimated  eten values to give levels 1 to 25 (from previous runs) :
etenf <- c(        21.40560, 25.32906, 28.86371, 32.07434, 35.01080, 37.71234 )
etenf <- c( etenf, 40.21043, 42.53071, 44.69434, 46.71898, 48.61956, 50.40878 )
etenf <- c( etenf, 52.08246, 53.35279, 54.56258, 55.83396, 57.16346, 58.54585 )
etenf <- c( etenf, 59.97358, 61.43734, 62.99543, 65.49976, 68.71259, 72.92473 )
etenf <- c( etenf, 78.62524  )  #  etenf values for Levels 1 to 25

 eten <- etenf
#  Unit vector of same length as eten
 eee  <-  1 + 0*eten        
#  ages at start of groups; later extended.
xx   <-  c(0, 1, seq(5,75, by=5) )

#   FEMALE ax and bx are the coefficients of the simple regressions of
#   the nqx's on e10
ax <- c(0.47504, 0.45025, 0.9376, 0.10041,0.10126, 0.11261, 0.13137, 0.15448, 0.17693, 0.18440, 0.19440, 0.22364, 0.30043, 0.41033, 0.56691, 0.77206, 0.96175)

bx <- c(-0.006923, -0.006805, -0.002928, -0.001497, -0.001480, -0.001618, -0.001893, -0.002239, -0.002566, -0.002612, -0.002712, -0.003011, -0.004053, -0.005394, -0.007187, -0.009334, -0.010681)


# a1x and b1x are the coefficients of the log-regressions of the nqx's on e10
a1x <- c(5.7332, 7.6298, 7.1271, 6.1089, 5.4984, 5.2649, 5.2547, 5.3691, 5.3186, 4.9099, 4.6164, 4.3673, 4.4363, 4.4163, 4.4030, 4.3826, 4.3108)

b1x <-  c(-0.05133, -0.08909, -0.08647, -0.07192, -0.05955, -0.5372, -0.05236, -0.05339, -0.05136,  -0.04261, -0.03627, -0.02961, -0.02858, -0.02511, -0.02152, -0.01784, -0.01355)

#####################################################################
##### #    Calculate an array of probabilities of dying  nqx  up to age 75
#######################################################################

#       The first subscript refers to e10; the second to the age x.
  nqxa <-    eee %o% ax  +  eten %o% bx                 #  nqx from linear model
  nqxl <-   (1/10000)*10^( eee %o% a1x   +  eten %o% b1x  )  # exponential model
#
#   Choose among   linear,  mixed, and exponential models
#       on basis of crossover points and slopes.

  slopechk <-  ( (0*eten + 1/log(10)) %o% (bx/b1x) ) *  (1/nqxl ) # Slope check
  nqx  <-  ( nqxa + nqxl) / 2
  nqx  <-  ifelse(  nqxa <= nqxl  &  slopechk <  1,  nqxa,  nqx )
  nqx  <-  ifelse(  nqxa <= nqxl  &  slopechk >= 1,  nqxl,  nqx )


####################################################################
####    Calculate the survivorship proportions   lx  up to age 75
##################################################################

  u   <-  log( 1 - nqx)
  uu  <-  t( apply( u, c(1), cumsum ) )
  km1 <-  length( uu[1,]) - 1    #  number of age groups minus 1
  lx  <-  cbind(  1+0*eee,  exp( uu[ , 1:km1])  )



#####################################################################
#   Special extension from ages 80 to 100 (i.e. for j> 18
###############################################################

     nL75 <-   5*lx[,17]*(1-nqx[,17])  + (1/2)*5*lx[,17]*nqx[,17]
     leighty <- lx[,17]*(1-nqx[,17])         #  Survivorship to age 80
     m77  <-   lx[,17]*nqx[,17]/nL75        #  m77 is hazard age 77.5.
     m105 <-   0.613+ 1.75*nqx[,17]
     kappa  <- (1/27.5)*log( m105 / m77 )  # THIS CORRECTS COALE'S MISTAKE
     m80    <-  m77*exp(kappa*2.5)           #  Hazard at age 80

     xlong <- seq(80,121, by = 0.2 ) - 80   #  ages in 1/5 years -80
     uu  <-    kappa %o% xlong
     u   <-    (-m80/kappa) %o% ( 1 + 0*xlong)
     llx  <-   leighty  %o% ( 1 + 0*xlong)
     llx  <-   llx  *  exp (  u*exp(uu) - u )

     nL80 <-  llx[,  1:26  ] %*% c( 1/2, rep(1,24), 1/2)*(1/5)
     nL85 <-  llx[, 26:51  ] %*% c( 1/2, rep(1,24), 1/2)*(1/5)
     nL90 <-  llx[, 51:76  ] %*% c( 1/2, rep(1,24), 1/2)*(1/5)
     nL95 <-  llx[, 76:101 ] %*% c( 1/2, rep(1,24), 1/2)*(1/5)
     T100 <-  as.vector( llx[, 101:201] %*% c( 1/2, rep(1,99), 1/2)*(1/5))


#################################################################
#####   Compute other functions in the lifetable
#################################################################


  xx   <-  c(0, 1, seq(5,95, by=5) )   #   ages
  nn   <-  c( diff(xx),8 )              #   age group widths
  lx   <-  cbind( lx[ , 1:17], llx[ , c(1,26,51,76)] )
  qqx  <-  ( llx[ , c(1,26,51) ] - llx[ , c(26,51,76)] )/llx[, c(1,26,51)]
  nqx  <-  cbind( nqx[ , 1:17], qqx, eee  )


  dimnames(nqx)[[1]] <- paste(seq(eten))
  dimnames(nqx)[[2]] <- paste(xx)
  dimnames(lx)       <- dimnames(nqx)

  na0  <-  0.050 + 3.00*nqx[,1]
  na0  <-  ifelse( nqx[,1]< 0.100, na0, 0.35 )
  na1  <-  1.524 - 1.625*nqx[,1]
  na1  <-  ifelse( nqx[,1]< 0.100, na1, 1.361)
  nax  <-  cbind( na0, na1, 2.25*eee)
  nax  <-  cbind( nax,  2.6*eee %o% rep(1,13), 2.5*eee)
  nax  <-  cbind( nax,  2.5*eee %o% rep(1,4)  )   #  Extend to length 21
  dimnames(nax)  <-  dimnames(nqx)

  ndx <-  lx*nqx
  nLx  <- ( eee %o% nn )*( lx - ndx ) +   nax*ndx
  nLx  <- cbind( nLx[ , 1:17], nL80, nL85, nL90, nL95 )
  nmx  <-  ndx/nLx
  K    <-  rev(seq(xx))         #  Reversed indices for age groups
  Tx   <-  t( apply(  nLx[, K],  c(1), cumsum) )
  Tx   <-  ( T100 %o%  rep(1,21) ) + Tx[,K]
  ex   <-  Tx/lx

  nLxf  <-  nLx        #   Save  nLx for females in separate array

  out <- list(age=xx,width=nn,lx=lx,nqx=nqx,nax=nax,ndx=ndx,nLx=nLx,nmx=nmx,Tx=Tx,ex=ex)

##################################################################
#####   Interpolate to find  eten values corresponding to Levels 1,2,3...
#####   This section is only used in first runs of program to determine eten
######################################################################

if(i==1){
   u  <-  approx(  ex[,1], eten,  seq(20, 80, by=2.5), rule=2 )
   etenf  <-  u$y      #  Estimate of eten which produce C-D LEvels
 }
  }


######################################################################


 if(sex=="F") return(out)
  else{
    #
    
# etenm is a linear function of etenf from runs of Coale-Demeny female
#       program etenm <- 2.298048 + 0.914067*etenf

 etenm <- c(       21.86420, 25.45051, 28.68141, 31.61614, 34.30026, 36.76965 )
 etenm <- c(etenm, 39.05308, 41.17397, 43.15167, 45.00232, 46.73959, 48.37505 )
 etenm <- c(etenm, 49.90490, 51.06607, 52.17190, 53.33403, 54.54928, 55.81287 )
 etenm <- c(etenm, 57.11792, 58.45590, 59.88009, 62.16922, 65.10596, 68.95614 )
 etenm <- c(etenm, 74.16678  )

## eten_seq(20,86, by = 2)   #  Code used for first pass to find etenm values
 eten <-  etenm               #
 eee <-  1 + 0*eten           #  Unit vector of same length as eten
 xx   <-  c(0, 1, seq(5,75, by=5) )   #  ages at start of groups; later extended.

# MALE ax and bx are the coefficients of the simple regressions of the
# nqx's on e10
 
ax <- c(0.54327, 0.46169, 0.18983, 0.09551, 0.09666, 0.13472, 0.14325, 0.15280, 0.17535, 0.20924, 0.24673, 0.28578, 0.36171, 0.45849, 0.59986, 0.82662, 1.03681) 

bx <-  c(-0.008251, -0.007290, -0.002974, -0.001476, -0.001422, -0.001968, -0.002103, -0.002244, -0.002589, -0.003083, -0.003605, -0.004016, -0.005037, 0.006124,-0.007677, -0.010241, -0.011906)

#  MALE a1x and b1x are the coefficients of the log-regressions of the
#  nqx's on e10

a1x <- c(5.6151, 7.2025, 6.1947, 5.3488, 4.5662, 4.6970, 4.7661, 4.7248, 4.7568, 4.7280, 4.6020, 4.3499, 4.3718, 4.2977, 4.2858, 4.3482, 4.3197)

b1x <-  c(-0.05022, -0.08475, -0.07195, -0.06047, -0.04322, -0.04277, -0.04372, -0.04236, -0.04197, -0.03986, -0.03578, -0.02857, -0.02682, -0.02244, -0.01913, -0.01710, -0.01357)

#####################################################################
##### #    Calculate an array of probabilities of dying  nqx  up to age 75
#######################################################################

#       The first subscript refers to e10; the second to the age x.
  nqxa <-    eee %o% ax  +  eten %o% bx                 #  nqx from linear model
  nqxl <-   (1/10000)*10^( eee %o% a1x   +  eten %o% b1x  )  # exponential model
#
#   Choose among   linear,  mixed, and exponential models
#       on basis of crossover points and slopes.

  slopechk <-  ( (0*eten + 1/log(10)) %o% (bx/b1x) ) *  (1/nqxl ) # Slope check
  nqx  <-  ( nqxa + nqxl) / 2
  nqx  <-  ifelse(  nqxa <= nqxl  &  slopechk <  1,  nqxa,  nqx )
  nqx  <-  ifelse(  nqxa <= nqxl  &  slopechk >= 1,  nqxl,  nqx )


####################################################################
####    Calculate the survivorship proportions   lx  up to age 75
##################################################################

  u   <-  log( 1 - nqx)
  uu  <-  t( apply( u, c(1), cumsum ) )
  km1 <-  length( uu[1,]) - 1    #  number of age groups minus 1
  lx  <-  cbind(  1+0*eee,  exp( uu[ , 1:km1])  )



#####################################################################
#   Special extension from ages 80 to 100 (i.e. for j> 18
###############################################################

     nL75 <-   5*lx[,17]*(1-nqx[,17])  + (1/2)*5*lx[,17]*nqx[,17]
     leighty <- lx[,17]*(1-nqx[,17])         #  Survivorship to age 80
     m77  <-   lx[,17]*nqx[,17]/nL75        #  m77 is hazard age 77.5.
     m105 <-   0.613+ 1.75*nqx[,17]
     kappa  <- (1/27.5)*log( m105 / m77 )  # THIS CORRECTS COALE'S MISTAKE
     m80    <-  m77*exp(kappa*2.5)           #  Hazard at age 80

     xlong <- seq(80,121, by = 0.2 ) - 80   #  ages in 1/5 years -80
     uu  <-    kappa %o% xlong
     u   <-    (-m80/kappa) %o% ( 1 + 0*xlong)
     llx  <-   leighty  %o% ( 1 + 0*xlong)
     llx  <-   llx  *  exp (  u*exp(uu) - u )

     nL80 <-  llx[,  1:26  ] %*% c( 1/2, rep(1,24), 1/2)*(1/5)
     nL85 <-  llx[, 26:51  ] %*% c( 1/2, rep(1,24), 1/2)*(1/5)
     nL90 <-  llx[, 51:76  ] %*% c( 1/2, rep(1,24), 1/2)*(1/5)
     nL95 <-  llx[, 76:101 ] %*% c( 1/2, rep(1,24), 1/2)*(1/5)
     T100 <-  as.vector( llx[, 101:201] %*% c( 1/2, rep(1,99), 1/2)*(1/5))


#################################################################
#####   Compute other functions in the lifetable
#################################################################


  xx   <-  c(0, 1, seq(5,95, by=5) )   #   ages
  nn   <-  c( diff(xx),8 )              #   age group widths
  lx   <-  cbind( lx[ , 1:17], llx[ , c(1,26,51,76)] )
  qqx  <-  ( llx[ , c(1,26,51) ] - llx[ , c(26,51,76)] )/llx[, c(1,26,51)]
  nqx  <-  cbind( nqx[ , 1:17], qqx, eee  )


  dimnames(nqx)[[1]] <- paste(seq(eten))
  dimnames(nqx)[[2]] <- paste(xx)
  dimnames(lx)       <- dimnames(nqx)

  na0  <-  0.050 + 3.00*nqx[,1]
  na0  <-  ifelse( nqx[,1]< 0.100, na0, 0.35 )
  na1  <-  1.524 - 1.625*nqx[,1]
  na1  <-  ifelse( nqx[,1]< 0.100, na1, 1.361)
  nax  <-  cbind( na0, na1, 2.25*eee)
  nax  <-  cbind( nax,  2.6*eee %o% rep(1,13), 2.5*eee)
  nax  <-  cbind( nax,  2.5*eee %o% rep(1,4)  )   #  Extend to length 21
  dimnames(nax)  <-  dimnames(nqx)

  ndx <-  lx*nqx
  nLx  <- ( eee %o% nn )*( lx - ndx ) +   nax*ndx
  nLx  <- cbind( nLx[ , 1:17], nL80, nL85, nL90, nL95 )
  nmx  <-  ndx/nLx
  K    <-  rev(seq(xx))         #  Reversed indices for age groups
  Tx   <-  t( apply(  nLx[, K],  c(1), cumsum) )
  Tx   <-  ( T100 %o%  rep(1,21) ) + Tx[,K]
  ex   <-  Tx/lx

  nLxm  <-  nLx        #   Save  nLx for males in separate array

  out <- list(age=xx,width=nn,lx=lx,nqx=nqx,nax=nax,ndx=ndx,nLx=nLx,nmx=nmx,Tx=Tx,ex=ex)

 return(out)
}

}

