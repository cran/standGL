ridgedProj <-
function(U,D,y){
  iprod <- rep(0,ncol(U))
  p <- rep(0,length(y))

  iprod <- D*(t(U) %*% y)
  p <- U %*% iprod
  return(p)
}

