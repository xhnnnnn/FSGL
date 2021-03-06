% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictFSGL.R
\name{predictFSGL}
\alias{predictFSGL}
\title{Outputs Predicted Responses from an FSGL Model for New Observations}
\usage{
predictFSGL(x, newX, lam)
}
\arguments{
\item{x}{fitted "FSGL" object}

\item{newX}{covariate matrix for new observations whose responses we wish to predict}

\item{lam}{the index of the lambda value for the model with which we desire to predict}
}
\description{
Outputs predicted response values for new user input observations at a specified \code{lambda} value
}
\details{
Predicted outcomes are given
}
\examples{
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
data = list(x = X, y = y)
Fit = FSGL(data, index, type = "linear")
X.new = matrix(rnorm(n * p), ncol = p, nrow = n)
predictFSGL(Fit, X.new, 5)
}
