logitProj <-
function(U,y){
  iprod <- rep(0,ncol(U))
  p <- rep(0,length(y))

  iprod <- t(U) %*% y
  p <- U %*% iprod
  return(p)
}

