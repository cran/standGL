\name{standGL}
\alias{standGL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits the standardized group lasso along the regularization path
}
\description{
Gives the solution to the standardized/ridged group lasso for linear and logistic regression for a number of values along the regularization path.
}
\usage{
standGL(y,X, index, family = "linear", alpha = 0.95, maxit = 1000, thresh = 10^(-3), nlam = 100, min.lam.frac = 0.05, is.pen = rep(1,length(index)), lam.path = NULL)

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
The sequence of models along the regularization path is fit by group descent
}
\value{
A list with $2$ elements
  \item{beta}{
a $p$ \code{x} \code{nlam} matrix of coefficient values, with column $i$ corresponding to the coefficient values for model fit with penalty parameter \code{lam.path}[$i$].
}
  \item{lam.path}{
the regularization path for which the model was fit (the actual sequence of \code{lambda} used).
}
}
\references{Simon, N. and Tibshirani, R. (2011)
  \emph{Standardization and the Group Lasso Penalty}, \cr
  \url{http://www-stat.stanford.edu/~nsimon/standGL.pdf}
  }
\author{Noah Simon and Rob Tibshirani\cr
Maintainer: Noah Simon \email{nsimon@stanford.edu}
}
\seealso{\code{cv.standGL}
}
\examples{
n = 100; p = 300; size.groups = 5
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 2*rnorm(n)
fit = standGL(y, X, index, family = "linear")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
