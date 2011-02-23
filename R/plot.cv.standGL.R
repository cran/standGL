plot.cv.standGL = function(x,...){
  cvobj = x
  xlab = "log(Lambda)"
  ylab = "CV Negative Log Likelihood"
  cvUp <- cvobj$cvNegLL + cvobj$cvNegLLsd
  cvDn <- cvobj$cvNegLL - cvobj$cvNegLLsd
  plot.args = list(x=log(cvobj$lambdas), y = cvobj$cvNegLL, xlab = xlab, ylab = ylab, ylim = range(cvUp, cvDn))
  do.call("plot", plot.args)
  error.bars(log(cvobj$lambdas), cvUp, cvDn, width = 0.01, col ="darkgrey")
  points(log(cvobj$lambdas),cvobj$cvNegLL,pch=20,col="red")
}
