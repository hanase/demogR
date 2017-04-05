`gen.time` <-
function(A,peryear=5){
  ro <- calc.ro(A)
  ea <- eigen.analysis(A)
  T <- peryear*log(ro)/log(ea$lambda)

  T
}

