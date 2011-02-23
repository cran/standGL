logitcalcLam <- function(y, Us, index, nlam = 20, min.lam = 0.05){
	groups <- unique(index)
	fitNorm <- rep(0,length(groups))

	for(i in groups){
		newInd <- which(index == i)
		fitNorm[i] <- sum((Us[[i]] %*% (t(Us[[i]]) %*% (y - mean(y)))^2))/length(newInd)
	}
	lam.max <- sqrt(max(fitNorm))
	lam.min <- min.lam * lam.max
	lam.list <- nlam:1
	scale <- (log(lam.max)-log(lam.min))/(nlam - 1)
	shift <- log(lam.min) - scale
	lam.list <- exp(lam.list * scale + shift)
	return(lam.list)
}

