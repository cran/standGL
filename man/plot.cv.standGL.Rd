\name{plot.cv.standGL}
\alias{plot.cv.standGL}
\title{plots the cross-validated error curve produced by cv.standGL}
\description{
Plots the cross-validated error curve, and confidence bounds for each \code{lambda} in our regularization path.}
\usage{
\method{plot}{cv.standGL}(x, ...)
}
\arguments{
  \item{x}{fitted \code{"cv.standGL"} object}
  \item{...}{additional arguments to be passed to plot}
}
\details{
A plot is produced}
\references{Simon, N. and Tibshirani, R. (2011)
  \emph{Standardization and the Group Lasso Penalty}, \cr
  \url{http://www-stat.stanford.edu/~nsimon/standGL.pdf}
  }
\author{Noah Simon and Rob Tibshirani\cr
Maintainer: Noah Simon \email{nsimon@stanford.edu}
}
\seealso{\code{standGL} and \code{cv.standGL}.}
\examples{
n = 100; p = 300; size.groups = 5
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 2*rnorm(n)
cvFit = cv.standGL(y, X, index, family = "linear")
plot(cvFit)
}
\keyword{models}
\keyword{regression}
