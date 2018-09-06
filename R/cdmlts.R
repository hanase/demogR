# Model South -------------------------------------
# With correction by Jim Oeppen (September 2018)
cdmlts <- function(sex = "F"){
    if (sex != "F" & sex !="M")
        stop("sex must be either F or M!")
  
    if(sex == "F") { # Female
        # CD lifetables indexed by e_10 = eten
        # Level 1 is ef=20,  Level 25 is  ef = 80
        # JO optimisation results
        etenf <- c(30.27658, 33.47231, 36.35843, 38.98633, 41.39559,
                   43.61651, 45.67469, 47.58958, 49.37782, 50.89556,
                   52.15870, 53.44483, 54.75041, 56.07245, 57.40630,
                   58.74726, 60.09140, 61.44296, 63.18308, 65.31128,
                   67.73256, 70.52334, 73.79105, 77.68029, 82.38960)

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
        # JO corrected; see Gough (1987) East and South reversed
        na1 <- 1.402 - 1.627 * nqx[, 1]
        na1 <- ifelse(nqx[, 1] < 0.1, na1, 1.324)

        nax  <-  cbind( na0, na1, 2.25*eee)
        nax  <-  cbind( nax,  2.6*eee %o% rep(1,13), 2.5*eee)
        nax  <-  cbind( nax,  2.5*eee %o% rep(1,4)  )   #  Extend to length 21
        dimnames(nax)  <-  dimnames(nqx)

        ndx <-  lx*nqx
        nLx  <- ( eee %o% nn )*( lx - ndx ) +   nax*ndx
        nLx  <- cbind( nLx[ , 1:17], nL80, nL85, nL90, nL95 )
        
        # JO add dimnames
        dimnames(nLx)[[1]] <- paste(seq(eten))
        dimnames(nLx)[[2]] <- paste(xx)
        
        nmx  <-  ndx/nLx
        K    <-  rev(seq(xx))         #  Reversed indices for age groups
        Tx   <-  t( apply(  nLx[, K],  c(1), cumsum) )
        Tx   <-  ( T100 %o%  rep(1,21) ) + Tx[,K]
        ex   <-  Tx/lx
        
        out <- list(age=xx,width=nn,lx=lx,nqx=nqx,nax=nax,ndx=ndx,nLx=nLx,nmx=nmx,Tx=Tx,ex=ex)

        return(out)
    } # end Female part

    # Male
    
    # JO optimisation results
    etenm <- c(31.56466, 34.27747, 36.72772, 38.95883, 41.00415,
               42.89005, 44.63729, 46.26285, 47.78132, 49.06971,
               50.14219, 51.23407, 52.34256, 53.46486, 54.59741,
               55.73606, 56.87751, 58.02539, 59.50311, 61.31097,
               63.36639, 65.73631, 68.51146, 71.81394, 75.81369)

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
     # JO corrected
     m105 <- 0.551 + 1.75 * nqx[, 17] # Male

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
    nn   <-  c( diff(xx),8 )             #   age group widths
    lx   <-  cbind( lx[ , 1:17], llx[ , c(1,26,51,76)] )
    qqx  <-  ( llx[ , c(1,26,51) ] - llx[ , c(26,51,76)] )/llx[, c(1,26,51)]
    nqx  <-  cbind( nqx[ , 1:17], qqx, eee  )

    dimnames(nqx)[[1]] <- paste(seq(eten))
    dimnames(nqx)[[2]] <- paste(xx)
    dimnames(lx)       <- dimnames(nqx)

    # JO corrected na0 and na1
    na0 <- 0.0425 + 2.875 * nqx[, 1]
    na0 <- ifelse(nqx[, 1] < 0.1, na0, 0.33)
    # JO corrected; see Gough (1987) East and South reversed
    na1 <- 1.541 - 3.013 * nqx[, 1]
    na1 <- ifelse(nqx[, 1] < 0.1, na1, 1.313)

    nax  <-  cbind( na0, na1, 2.25*eee)
    nax  <-  cbind( nax,  2.6*eee %o% rep(1,13), 2.5*eee)
    nax  <-  cbind( nax,  2.5*eee %o% rep(1,4)  )   #  Extend to length 21
    dimnames(nax)  <-  dimnames(nqx)

    ndx <-  lx*nqx
    nLx  <- ( eee %o% nn )*( lx - ndx ) +   nax*ndx
    nLx  <- cbind( nLx[ , 1:17], nL80, nL85, nL90, nL95 )
    # JO add dimnames
    dimnames(nLx)[[1]] <- paste(seq(eten))
    dimnames(nLx)[[2]] <- paste(xx)
    
    nmx  <-  ndx/nLx
    K    <-  rev(seq(xx))         #  Reversed indices for age groups
    Tx   <-  t( apply(  nLx[, K],  c(1), cumsum) )
    Tx   <-  ( T100 %o%  rep(1,21) ) + Tx[,K]
    ex   <-  Tx/lx

    out <- list(age=xx,width=nn,lx=lx,nqx=nqx,nax=nax,ndx=ndx,nLx=nLx,nmx=nmx,Tx=Tx,ex=ex)

    return(out)
}

