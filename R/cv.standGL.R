cv.standGL <- function(y,X, index, family = "linear", alpha = 0.95, maxit = 1000, thresh = 10^(-3), nlam = 100, nfold = 10, min.lam.frac = 0.05, is.pen = rep(1,length(index)), lam.path = NULL){

  if(family == "linear"){
    output <- cv.ridgedGL(y, X, index, lam.path = lam.path, thresh = thresh, maxit = maxit, nfold = nfold, nlam = nlam, min.lam.frac = min.lam.frac, alpha = 1, is.pen = is.pen)
}

  if(family == "ridge"){
    output <- cv.ridgedGL(y, X, index, lam.path = lam.path, thresh = thresh, maxit = maxit, nfold = nfold, nlam = nlam, min.lam.frac = min.lam.frac, alpha = alpha, is.pen = is.pen)
}

  if(family == "logit"){
    output <- cv.logitGL(y, X, index, lam.path = lam.path, thresh = thresh, maxit = maxit, nlam = nlam, min.lam.frac = min.lam.frac, is.pen = is.pen)
}

return(output)
}
