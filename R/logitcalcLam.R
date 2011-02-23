logitcalcLam <- function(y, Us, index, nlam = 20, min.lam = 0.05, is.pen = rep(1,length(unique(index)))){
  groups <- unique(index)
  fitNorm <- rep(0,length(groups))
  
  for(i in 1:length(groups)){
    if(is.pen[i] != 0){
      newInd <- which(index == groups[i])
      fitNorm[i] <- sum((Us[[i]] %*% (t(Us[[i]]) %*% (y - 0.5)))^2)/length(newInd)
    }
  }
  lam.max <- sqrt(max(fitNorm))
  lam.min <- min.lam * lam.max
  lam.list <- nlam:1
  scale <- (log(lam.max)-log(lam.min))/(nlam - 1)
  shift <- log(lam.min) - scale
  lam.list <- exp(lam.list * scale + shift)
  return(lam.list)
}

