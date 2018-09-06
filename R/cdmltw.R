# Model West ------------------------------------
# With correction by Jim Oeppen (September 2018)

cdmltw <- function(sex = "F"){
    if (sex != "F" & sex !="M")
        stop("sex must be either F or M!")
  
    # CD lifetables indexed by e_10 = eten

    if(sex == "F") { # Female
    # JO optimisation results
        etenf <- c(21.40549, 25.32904, 28.86351, 32.07396, 35.01032,
                   37.71173, 40.20992, 42.53003, 44.69388, 46.71830,
                   48.61905, 50.40814, 52.08187, 53.35198, 54.56192,
                   55.83317, 57.16274, 58.54485, 59.97243, 61.43625,
                   62.99381, 65.49800, 68.71056, 72.92253, 78.62273)
        eten <- etenf
        #  Unit vector of same length as eten
        eee  <-  1 + 0*eten        
        #  ages at start of groups; later extended.
        xx   <-  c(0, 1, seq(5,75, by=5) )
        
        #   FEMALE ax and bx are the coefficients of the simple regressions of
        #   the nqx's on e10
        ax <- c( .53774,.39368,.10927,.08548,.10979,.1358,.15134,.17032,.18464,.1939,.20138)
        ax <- c(ax,.2535,.31002,.43445,.53481,.69394,.84589)
        bx <- c(   -.008044,-.006162,-.001686,-.00132,-.001672,-.002051,-.002276,-.002556)
        bx <- c(bx,-.002745,-.002828,-.002831,-.003487,-.004118,-.005646,-.00646,-.007713)
        bx <- c(bx,-.008239)
        
        #     a1x and b1x are the coefficients of the log-regressions of the nqx's on e10.
        a1x <- c(    5.8992,7.4576,6.2018,5.9627,5.9335,5.9271,5.8145,5.6578,5.3632,4.96)
        a1x <- c(a1x,4.5275,4.4244,4.3131,4.3439,4.2229,4.1838,4.1294)
        b1x <- c(    -.05406,-.08834,-.0741,-.07181,-.06812,-.06577,-.06262,-.05875,-.05232)
        b1x <- c(b1x,-.0438,-.03436,-.03004,-.02554,-.02295,-.01773,-.01376,-.00978)
        
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
        # JO corrected originally 1.625; see Gough (1987)
        na1  <-  1.524 - 1.627*nqx[,1]
        na1  <-  ifelse( nqx[,1]< 0.100, na1, 1.361)
        nax  <-  cbind( na0, na1, 2.25*eee)
        nax  <-  cbind( nax,  2.6*eee %o% rep(1,13), 2.5*eee)
        nax  <-  cbind( nax,  2.5*eee %o% rep(1,4)  )   #  Extend to length 21
        dimnames(nax)  <-  dimnames(nqx)
        
        ndx <-  lx*nqx
        nLx  <- ( eee %o% nn )*( lx - ndx ) +   nax*ndx
        nLx  <- cbind( nLx[ , 1:17], nL80, nL85, nL90, nL95 )
        # JO added dimnames
        dimnames(nLx)[[1]] <- paste(seq(eten))
        dimnames(nLx)[[2]] <- paste(xx)
        
        nmx  <-  ndx/nLx
        K    <-  rev(seq(xx))         #  Reversed indices for age groups
        Tx   <-  t( apply(  nLx[, K],  c(1), cumsum) )
        Tx   <-  ( T100 %o%  rep(1,21) ) + Tx[,K]
        ex   <-  Tx/lx
        
        out <- list(age = xx, width = nn, e10 = eten, 
                    lx=lx,nqx=nqx,nax=nax,ndx=ndx,nLx=nLx,nmx=nmx,Tx=Tx,ex=ex)
        
        return(out)
    } # end Female
    
    # Male
    
    # JO optimisation results
    etenm <- c(21.86761, 25.44740, 28.67272, 31.60223, 34.28159,
               36.74631, 39.02583, 41.14277, 43.11689, 44.96458,
               46.69839, 48.33115, 49.85796, 51.01729, 52.12115,
               53.28135, 54.49441, 55.75599, 57.05851, 58.39447,
               59.81638, 62.10227, 65.03413, 68.87841, 74.08120)

    eten <-  etenm               #
    eee <-  1 + 0*eten           #  Unit vector of same length as eten
    xx   <-  c(0, 1, seq(5,75, by=5) )   #  ages at start of groups; later extended.

    #   MALE ax and bx are the coefficients of the simple regressions of the nqx's on e10.
    ax <- c(  0.63726, 0.40548, 0.10393, 0.07435, 0.09880, 0.14009, 0.15785, 0.18260, 0.21175)
    ax <- c( ax, 0.25049, 0.27894, 0.33729, 0.38425, 0.48968, 0.59565, 0.73085, 0.98976)
    bx <- c(    -0.009958, -0.006653, -0.001662, -0.001183, -0.001539, -0.002183, -0.002479)
    bx <- c(bx, -0.002875, -0.003312, -0.003864, -0.004158, -0.004856, -0.005190, -0.006300)
    bx <- c(bx, -0.007101, -0.007911, -0.008695 )

    #    MALE  a1x and b1x are the coefficients of the log-regressions of the nqx's on e10.
    a1x <- c(      5.8061, 7.1062, 5.4472, 5.0654, 4.8700, 5.0677, 5.2660, 5.3438, 5.2792)
    a1x <- c(a1x,  5.0415, 4.6666, 4.4506, 4.2202, 4.1851, 4.1249, 4.1051, 4.1133)
    b1x <- c(     -0.05338, -0.08559, -0.06295, -0.05817, -0.05070, -0.05156, -0.05471)
    b1x <- c(b1x, -0.05511, -0.05229, -0.04573, -0.03637, -0.02961, -0.02256 )
    b1x <- c(b1x, -0.01891, -0.01491, -0.01161, -0.00895 )

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
     # JO corrected (0.551 instead of 0.613)
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
    nn   <-  c( diff(xx),8 )              #   age group widths
    lx   <-  cbind( lx[ , 1:17], llx[ , c(1,26,51,76)] )
    qqx  <-  ( llx[ , c(1,26,51) ] - llx[ , c(26,51,76)] )/llx[, c(1,26,51)]
    nqx  <-  cbind( nqx[ , 1:17], qqx, eee  )

    dimnames(nqx)[[1]] <- paste(seq(eten))
    dimnames(nqx)[[2]] <- paste(xx)
    dimnames(lx)       <- dimnames(nqx)
    
    # JO corrected na0 and na1
    na0 <- 0.0425 + 2.875 * nqx[, 1]
    na0 <- ifelse(nqx[, 1] < 0.1, na0, 0.33)
    na1 <- 1.653 - 3.013 * nqx[, 1]
    na1 <- ifelse(nqx[, 1] < 0.1, na1, 1.352)

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

    out <- list(age = xx,width = nn, e10 = eten, 
                lx=lx,nqx=nqx,nax=nax,ndx=ndx,nLx=nLx,nmx=nmx,Tx=Tx,ex=ex)

    return(out)
}

