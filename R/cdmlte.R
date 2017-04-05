# modified version of Ken Wachter's routine

`cdmlte` <-
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
ax <- c(0.78219, 0.46584, 0.13739, 0.07600, 0.10067, 0.13039, 0.15401, 0.16941, 0.18184, 0.18555, 0.19407, 0.24415, 0.34490, 0.49585, 0.68867, 0.88452, 1.07727)

bx <-  c(-0.011679,  -0.007284, -0.002136, -0.001166, -0.001529, -0.001973, -0.002335, -0.002559, -0.002718, -0.002718, -0.002746, -0.003376, -0.004723, -0.006651, -0.008874, -0.010551, -0.011513)

# a1x and b1x are the coefficients of the log-regressions of the nqx's on e10
a1x <- c(5.8529, 7.2269, 6.3204, 5.6332, 5.5780, 5.5872, 5.6149, 5.4593, 5.1881, 4.8186, 4.4509, 4.3702, 4.4480, 4.4917, 4.4702, 4.3759, 4.2972)

b1x <-  c(-0.05064, -0.08351, -0.07590, -0.06684, -0.06295, -0.06081, -0.06004, -0.05616, -0.05000, -0.04209, -0.03368, -0.02966, -0.02807, -0.02544, -0.02152, -0.01640, -0.01191)

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
 
ax <- c(1.07554, 0.55179, 0.15292, 0.06856, 0.10060, 0.14725, 0.15127, 0.17022, 0.20786, 0.24876, 0.28685, 0.32623, 0.38906, 0.49337, 0.66168, 0.84188, 1.03876)

bx <-  c(-0.017228, -0.009201, -0.002523, -0.001096, -0.001578, -0.002312, -0.002381, -0.002686, -0.003277,-0.003868, -0.004320, -0.004654, -0.005243, -0.006341, -0.008182, -0.009644, -0.010780)

#  MALE a1x and b1x are the coefficients of the log-regressions of the
#  nqx's on e10

a1x <- c(6.3796, 7.8944, 6.4371, 5.1199, 4.9229, 5.1056, 5.1036, 5.1685, 5.1986, 5.0221, 4.6915, 4.3492, 4.1849, 4.1647, 4.2175, 4.2171, 4.2155)

b1x <-  c(-0.06124, -0.09934, -0.0876, -0.05978, -0.05182, -0.05225, -0.05207, -0.05244, -0.05131, -0.04577, -0.03697, -0.02767, -0.02171, -0.01842, -0.01634, -0.01324, -0.01035)
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


