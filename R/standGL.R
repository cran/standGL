standGL <- function(y,X, index, family = "linear", alpha = 0.95, maxit = 1000, thresh = 10^(-3), nlam = 100, min.lam.frac = 0.05, is.pen = rep(1,length(unique(index))), lam.path = NULL){

  if(family == "linear"){
    output <- ridgedGL(y, X, index, lam.path = lam.path, thresh = thresh, maxit = maxit, nlam = nlam, min.lam.frac = min.lam.frac, alpha = 1, is.pen = is.pen)
}

  if(family == "ridge"){
    output <- ridgedGL(y, X, index, lam.path = lam.path, thresh = thresh, maxit = maxit, nlam = nlam, min.lam.frac = min.lam.frac, alpha = alpha, is.pen = is.pen)
}

  if(family == "logit"){
    output <- logitGL(y, X, index, lam.path = lam.path, thresh = thresh, maxit = maxit, nlam = nlam, min.lam.frac = min.lam.frac, is.pen = is.pen)
}

return(output)
}
