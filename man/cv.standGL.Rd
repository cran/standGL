\name{cv.standGL}
\alias{cv.standGL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finds the optimal penalty parameter along the regularization path for the standardized group lasso
}
\description{
Cross Validates to find the optimal penalty parameter, lambda, for the standardized (and ridged) group lasso for linear and logistic regression.
}
\usage{
cv.standGL(y,X, index, family = "linear", alpha = 0.95, maxit = 1000, thresh = 10^(-3), nlam = 100, nfold = 10, min.lam.frac = 0.05, is.pen = rep(1,length(unique(index))), lam.path = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response variable, should be an n-vector (where n is the number of obs), for logistic regression this is a vector of 1s and 0s.
}
  \item{X}{
input matrix, n by p (where n is the number of obs, and p is the number of covariates).
}
  \item{index}{
a p-vector of group memberships for the variables (used in the penalty)
}

  \item{family}{
model type: one of ("ridge","linear","logit").  "ridge" corresponds to a ridged linear model (a ridged logit function has not been coded in)
}

  \item{alpha}{
The ridged mixing parameter. \code{alpha} = 1 is the standardized group lasso penalty.  This is only relevant for the "ridged" family (in the literature, \code{delta} is used, here \code{delta} = (1-\code{alpha})*\code{lambda}).
}

  \item{maxit}{
the maximum number of passes over all the groups, default is 100
}

  \item{thresh}{
the convergence threshhold for groupwise descent, default is 10^(-4)
}

  \item{nlam}{
number of lambdas in the regularization path, default is 100
}

  \item{nfold}{
number of folds in the cross validation loop, default is 10.
}

  \item{min.lam.frac}{
smallest value of the minimum lambda, as a function of maximum lambda, default is 0.05.  This value is irrelevant if a user specificed lam.path is given.
}

  \item{is.pen}{
a 0-1 vector indicating whether or not each group is to be penalized (one should use an unpenalized column of 1s, coded as its own group as the intercept)
}

  \item{lam.path}{
an optional user specified regularization path.  If left empty, the package will calculate its own path
}

}
\details{
The function runs \code{standGL} \code{nfold}+1 times; the initial run is to find the \code{lambda} sequence, subsequent runs are used to compute the cross-validated error rate and its standard deviation.
}
\value{
An object with S3 class \code{"cv.standGL"}
  \item{cvNegLL}{
An \code{nlam} vector of cross validated negative log likelihoods (squared error loss in the \code{linear} case, along the regularization path)
}
  \item{cvNegLLsd}{
An \code{nlame} vector of approximate standard deviations of \code{cdNegLL}
}
  \item{lambdas}{The actual list of \code{lambda} values used in the regularization path.
}
  \item{lambda.min}{The \code{lambda} value which achieves the minimum cross validated error.
}
  \item{fit}{A model fit object created by a call to \code{standGL} on the entire dataset
}
}
\references{Simon, N. and Tibshirani, R. (2011)
  \emph{Standardization and the Group Lasso Penalty}, \cr
  \url{http://www-stat.stanford.edu/~nsimon/standGL.pdf}
  }
\author{Noah Simon and Rob Tibshirani\cr
Maintainer: Noah Simon \email{nsimon@stanford.edu}
}
\seealso{\code{standGL}
}
\examples{
set.seed(1)
n = 30; p = 50; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 2*rnorm(n)
cvFit = cv.standGL(y, X, index, family = "linear")
}
\keyword{models}
\keyword{regression}
