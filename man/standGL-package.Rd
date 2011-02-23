\name{standGL-package}
\alias{standGL-package}
\docType{package}
\title{
Fits the standardized (and ridged) group lasso penalty for linear and logistic regression
}
\description{
This packages fits and cross validates the standardized (and ridged) group lasso for linear and logistic regression via groupwise descent.
}
\details{
\tabular{ll}{
Package: \tab standGL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-02-23\cr
License: \tab GNU\cr
LazyLoad: \tab yes\cr
}
Very straightforward to use. Data is input in the form \code{y,X}, and produces solutions over a path of \code{lambda} values.  Three functions:
\code{standGL} \cr
\code{cv.standGL} \cr
\code{plot.cv.standGL}
}
\author{
Noah Simon and Rob Tibshirani \cr

Maintainer: Noah Simon <nsimon@stanford.edu>
}
\references{
Simon, N. and Tibshirani, R. (2011)
  \emph{Standardization and the Group Lasso Penalty}, \cr
  \url{http://www-stat.stanford.edu/~nsimon/standGL.pdf}
}
\keyword{models}
\keyword{regression}
\keyword{ package }
\examples{
n = 30; p = 50; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 2*rnorm(n)
Fit = standGL(y, X, index, family = "linear")
cvFit = cv.standGL(y, X, index, family = "linear")
plot(cvFit)
}
