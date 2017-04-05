`coale` <-
function(b1,b4,nMx){
   if(nMx[1]>0.107){
     b1 <- c(0.350,0)
     b4 <- c(1.361,0)
   }
   nax12 <- c(0,0)
   nax12[1] <- b1[1] + b1[2] *nMx[1]
   nax12[2] <- b4[1] + b4[2]* nMx[1]
   return(nax12)
 }

