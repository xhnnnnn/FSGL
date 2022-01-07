#' @title Outputs Predicted Responses from an FSGL Model for New Observations
#' @description Outputs predicted response values for new user input observations at a specified \code{lambda} value
#' @param x fitted "FSGL" object
#' @param newX covariate matrix for new observations whose responses we wish to predict
#' @param lam the index of the lambda value for the model with which we desire to predict
#'
#' @details Predicted outcomes are given
#'
#' @export
#'
#' @examples
#' n = 50; p = 100; size.groups = 10
#' index <- ceiling(1:p / size.groups)
#' X = matrix(rnorm(n * p), ncol = p, nrow = n)
#' beta = (-2:2)
#' y = X[,1:5] %*% beta + 0.1*rnorm(n)
#' data = list(x = X, y = y)
#' Fit = FSGL(data, index, type = "linear")
#' X.new = matrix(rnorm(n * p), ncol = p, nrow = n)
#' predictFSGL(Fit, X.new, 5)
predictFSGL = function(x,newX,lam){
  cvobj = x

  X <- newX

  if(is.matrix(X)){
    X <- t(t(newX) - x$X.transform$X.means)
    if(!is.null(x$X.transform$X.scale)){
      X <- t(t(X) / x$X.transform$X.scale)
    }
  }
  if(is.vector(X)){
    X <- X - x$X.transform$X.means
    if(!is.null(x$X.transform$X.scale)){
      X <- X / x$X.transform$X.scale
    }
  }

  intercept <- x$intercept

  if(is.matrix(X)){
    eta <- X %*% x$beta[,lam] + intercept
  }
  if(is.vector(X)){
    eta <- sum(X * x$beta[,lam]) + intercept
  }

  y.pred <- eta

  return(y.pred)
}
