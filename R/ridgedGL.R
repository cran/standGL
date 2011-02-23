ridgedGL <-
function(y, X, index, lam.path = NULL, thresh = 10^(-4), maxit = 100, nlam = 100, min.lam.frac = 0.05, alpha = 0.5, is.pen = rep(1,length(index))){

  if(is.null(lam.path)){
    lam.path <- ridgedcalcLam(y, X, index, nlam, min.lam.frac)
  }

  numGroups <- length(unique(index))
  groupLen <- rep(0,numGroups)
  for(i in 1:numGroups){
    groupLen[i] <- length(which(index == unique(index)[i]))
  }
  
  Us <- list()
  Vs <- list()
  Ds <- list()

  for(i in 1:numGroups)
    {
      groupiMat <- X[,which(index == unique(index)[i])]
      svdDecomp <- svd(groupiMat)
      Us[[i]] <- svdDecomp$u
      Vs[[i]] <- svdDecomp$v
      Ds[[i]] <- svdDecomp$d
    }

  oldFit <- matrix(0, ncol = length(index), nrow = length(y))
  resid <- y

  Betas <- matrix(0, ncol = length(lam.path), nrow = ncol(X))
  
  isActive = rep(0, numGroups)
  
  for(i in 1:length(lam.path)){
    
    fit <- ridgedSingleFit(y, X, index, lam.path[i], thresh = thresh, maxit = maxit, Us = Us, Vs = Vs, Ds = Ds, alpha = alpha, oldFit = oldFit, resid = resid, isActive = isActive, is.pen = is.pen)

    resid <- fit$resid
    oldFit <- fit$oldFit
    isActive <- fit$isActive
    Betas[,i] <- fit$Betas
  }
  return(list(beta = Betas, lam.path = lam.path))
}
