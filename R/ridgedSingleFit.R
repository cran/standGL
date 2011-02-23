ridgedSingleFit <-
function(y, X, index, lambda = 1, thresh = 0.0001, maxit = 1000, Us, Vs, Ds, alpha, oldFit, resid, isActive, is.pen){

  numGroups <- length(unique(index))
  groupLen <- rep(0,numGroups)
  for(i in 1:numGroups){
    groupLen[i] <- length(which(index == unique(index)[i]))
  }

  newD <- list(c(1,2,3))  # stupid workaround...
  dfs <- rep(0,numGroups)

  for(j in 1:numGroups){
    newD[[j+1]]<- Ds[[j]]^2/(Ds[[j]]^2 + (1-alpha)*lambda*is.pen[j])
    dfs[j] <- sum(newD[[j+1]])
  }
  for(j in 1:numGroups){  # stupid workaround...
   newD[[j]] <- newD[[j+1]]
  }


  diff <- 1
  iter <- 1
  innerConverged <- 1
  useGroups <- 1:numGroups
  groupChange <- 1
  qqq <- 0

  while(groupChange == 1){

    groupChange <- 0
    innerConverged <- 1
    diff <- 1
    
    while(diff > thresh && iter < maxit){
      iter <- iter + 1
      diff <- 0

      if(innerConverged == 1){
        useGroups <- 1:numGroups
      }
      if(innerConverged == 0){
        useGroups <- which(isActive == 1)
      }
      
      for(i in useGroups){
        resid <- resid + oldFit[,i]
        proj <- ridgedProj(Us[[i]], newD[[i]], resid)
        normProj <- sqrt(sum((proj)^2))
        shrinkage <-  max(c((1-sqrt(dfs[i])*lambda*alpha*is.pen[i] / normProj),0))
        newFit <- shrinkage * proj
        resid <- resid - newFit
        diff <- diff + sum(abs(oldFit[,i] - newFit))
        
        if(max(abs(oldFit[,i])) == 0 && max(abs(newFit)) > 0){
          isActive[i] <- 1
          groupChange <- 1
        }
        oldFit[,i] <- newFit
      }
      innerConverged <- 0
      diff <- diff / nrow(X)
    }
  }

Betas <- list()
  
  for(i in 1:numGroups)
    {
      fit.ind <- which(abs(Ds[[i]]) > 10^(-5))
      beta <- Vs[[i]][,fit.ind] %*% ((1/Ds[[i]][fit.ind]) * t(Us[[i]]) %*% oldFit[,i]) 
      Betas[[i]] <- beta
    }
  Betas <- unlist(Betas)
  return(list(Betas = Betas, oldFit = oldFit, resid = resid, isActive = isActive))
}

