ridgedcalcLam <-
function(y, X, index, nlam = 20, min.lam = 0.05){
	groups <- unique(index)
	fitNorm <- rep(0,length(groups))
	normy <- sqrt(sum(y^2))

	for(i in groups){
		newInd <- which(index == i)
		if(length(newInd) >= nrow(X)){
                  fitNorm[i] <- normy
                }
		else{
		fitNorm[i] <- sqrt(sum((X[,newInd]%*%solve(t(X[,newInd])%*%X[,newInd])%*%t(X[,newInd]) %*% y)^2)/length(newInd))
		}
	}
	lam.max <- max(fitNorm)
	lam.min <- min.lam * lam.max
	lam.list <- nlam:1
	scale <- (log(lam.max)-log(lam.min))/(nlam - 1)
	shift <- log(lam.min) - scale
	lam.list <- exp(lam.list * scale + shift)
	return(lam.list)
}

