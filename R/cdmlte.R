# modified version of Ken Wachter's routine
# With correction by Jim Oeppen (September 2018)
cdmlte <- function(sex = "F"){
  if (sex != "F" & sex !="M")
    stop("sex must be either F or M!")
  
    if(sex == "F") { # Female
        # CD lifetables indexed by e_10 = eten
        # Level 1 is ef=20,  Level 25 is  ef = 80
        # JO optimisation output
        etenf <- c( 30.30972, 33.26908, 35.96734, 38.44387, 40.73029,
                    42.85097, 44.82695, 46.67433, 48.40755, 50.03848,
                    51.39305, 52.45719, 53.56302, 54.70946, 55.89505,
                    57.11504, 58.36643, 59.64167, 60.93385, 62.23599,
                    63.98995, 66.48732, 69.56091, 73.44775, 78.48869)

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

        # JO corrected na0 and na1
        na0 <- 0.01 + 3 * nqx[, 1]
        na0 <- ifelse(nqx[, 1] < 0.1, na0, 0.31)
        # JO corrected; see Gough (1987) East and South reversed
        na1 <- 1.487 - 1.627 * nqx[, 1]
        na1 <- ifelse(nqx[, 1] < 0.1, na1, 1.239)

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
    } # end Female

    # Male
    
    # JO optimisation output
    etenm <- c(33.09021, 35.44552, 37.59309, 39.56416, 41.38376,
               43.07166, 44.64422, 46.11452, 47.49416, 48.79216,
               49.87003, 50.71672, 51.59681, 52.50972, 53.45322,
               54.42448, 55.42053, 56.43549, 57.46425, 58.50057,
               59.89758, 61.88540, 64.33229, 67.42635, 71.43957)

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
    nn   <-  c( diff(xx),8 )              #   age group widths
    lx   <-  cbind( lx[ , 1:17], llx[ , c(1,26,51,76)] )
    qqx  <-  ( llx[ , c(1,26,51) ] - llx[ , c(26,51,76)] )/llx[, c(1,26,51)]
    nqx  <-  cbind( nqx[ , 1:17], qqx, eee  )

    dimnames(nqx)[[1]] <- paste(seq(eten))
    dimnames(nqx)[[2]] <- paste(xx)
    dimnames(lx)       <- dimnames(nqx)

    # JO corrected na0 and na1
    na0 <- 0.0025 + 2.875 * nqx[, 1]
    na0 <- ifelse(nqx[, 1] < 0.1, na0, 0.29)
    # JO corrected; see Gough (1987) East and South reversed
    na1 <- 1.614 - 3.013 * nqx[, 1]
    na1 <- ifelse(nqx[, 1] < 0.1, na1, 1.240)

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


