`cdmlts` <-
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
ax <- c(.52069,.68268,.17066,.09,.12189,.15083,.16073,.16719,.17408,.17278,.17800,.22639,.30167,.47682,.6744,.92943,1.16023)

bx <- c(-.007051,-.010453,-.002657,-.00138, -.001851,-.002279,-.002412,-.002505,-.002583,-.002504,-.002513,-.003140,-.004130,-.006501,-.008891,-.011532, -.013009)


# a1x and b1x are the coefficients of the log-regressions of the nqx's on e10
a1x <- c(4.5097,5.9815,5.6479,5.1045,5.2384,5.1708,5.0949,4.9291,4.8035,4.4917, 4.2693,4.1982,4.2724,4.4242,4.4554,4.4348,4.3542)

b1x <- c(-.02566,-.05332,-.06136,-.05537,-.05494,-.05171,-.04945,-.04590,-.0428,-.03615,-.03092,-.02717,-.02588,-.02491,-.02190,-.01775,-.01296)

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
 
ax <- c(.61903,.70613,.16455,.07634,.11449,.17104,.17171,.16483,.17905,.20606,.23208,.28,.35245,.49465,.66947,.89759,1.10111)

bx <- c( -.008974,-.011375,-.002674,-.001207,-.00181,-.002693,-.00271,-.002535,  -.002734,-.003081,-.003370,-.003917,-.004765,-.006569,-.008608, -.010843,-.011806)

#  MALE a1x and b1x are the coefficients of the log-regressions of the
#  nqx's on e10

a1x <- c(4.7096,6.3246,5.64,4.6816,4.9454,5.2748,5.1168,4.8459,4.766,4.5796,4.3559,4.1918,4.1492,4.2479,4.3069,4.3251,4.2684)

b1x <- c(-.0298,-.06433,-.06389,-.05008,-.05170,-.05458,-.05152,-.04547,-.04292,-.03738,-.03116,-.02547,-.02193,-.02063,-.01863,-.01552,-.01123 )


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

