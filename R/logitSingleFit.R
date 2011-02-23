logitSingleFit <-
function(y, X, index, lambda = 1, thresh = 0.0001, maxit = 1000, Us, Vs, oldFit, resid, isActive, is.pen){

  numGroups <- length(unique(index))
  groupLen <- rep(0,numGroups)
  for(i in 1:numGroups){
    groupLen[i] <- length(which(index == unique(index)[i]))
  }

  newD <- NULL
  dfs <- rep(0,numGroups) 

  outerIter <- 1
  outerDiff <- 1
  oldOuterFit <- apply(oldFit, 1, sum)

  while(outerDiff > thresh && outerIter < maxit){
    ex <- exp(oldOuterFit)
    probs <- ex/(1+ex)
#    resp <- oldOuterFit + 4 * (y-nrow(y)*probs)/(probs*(1-probs))
    resp <- oldOuterFit + 2 * (y-probs)
    resid <- 1*(y - probs)
    
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
          proj <- logitProj(Us[[i]], resid)
          normProj <- sqrt(sum((proj)^2))
          if(normProj != 0){
            shrinkage <-  max(c((1-sqrt(groupLen[i])*1*lambda*is.pen[i] / normProj),0)) ## 2 -> 4?
          }
          if(normProj == 0){
            shrinkage <- 0
          }
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
    newOuterFit <- apply(oldFit,1,sum)
    outerDiff <- mean(abs(newOuterFit - oldOuterFit))
    oldOuterFit <- newOuterFit
  }

Betas <- list()
  
  for(i in 1:numGroups)
    {
      beta <- Vs[[i]] %*% (t(Us[[i]]) %*% oldFit[,i]) 
      Betas[[i]] <- beta
    }
  Betas <- unlist(Betas)
  return(list(Betas = Betas, oldFit = oldFit, resid = resid, isActive = isActive))
}

