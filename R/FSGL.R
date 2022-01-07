#' @title Fit a Linear Model with a Combination of Lasso and Group Lasso Regularization
#' @description Fit a regularized generalized linear model via penalized maximum likelihood. The model is fit for a path of values of the penalty parameter. Fits linear, logistic and Cox models.
#'
#' @param data data should be a list of $x$ and $y$, x is a data matrix (n x p) and y is a vector
#' @param index A p-vector indicating group membership of each covariate
#' @param maxit Maximum number of iterations to convergence
#' @param thresh Convergence threshold for change in beta
#' @param min.frac The minimum value of the penalty parameter, as a fraction of the maximum value
#' @param nlam Number of lambda to use in the regularization path
#' @param gamma Fitting parameter used for tuning backtracking (between 0 and 1)
#' @param standardize Logical flag for variable standardization prior to fitting the model
#' @param verbose Logical flag for whether or not step number will be output
#' @param step Fitting parameter used for inital backtracking step size (between 0 and 1)
#' @param reset Fitting parameter used for taking advantage of local strong convexity in nesterov momentum (number of iterations before momentum term is reset)
#' @param alpha The mixing parameter. \code{alpha} = 1 is the lasso penalty. \code{alpha} = 0 is the group lasso penalty.
#' @param lambdas A user specified sequence of lambda values for fitting. We recommend leaving this NULL and letting SGL self-select values
#' @param UseUpperBound A logical flag for using upper bound
#'
#' @return
#' An object of type "FSGL"
#' \describe{
#'     \item{beta}{A p by $nlam$ matrix, giving the penalized MLEs for the nlam different models, where the index corresponds to the penalty parameter $lambda$}
#'     \item{lambdas}{The actual sequence of $lambda$ values used (penalty parameter)}
#'     \item{X.transform}{A list used in $predict$ which gives the empirical mean and variance of the x matrix used to build the model}
#' }
#' @export
#'
#' @examples
#' library(FSGL)
#'
#' n = 100; p = 20; size.groups = 10
#' index <- ceiling(1:p / size.groups)
#' X = matrix(rnorm(n * p), ncol = p, nrow = n)
#' beta = (-2:2)
#' y = X[,1:5] %*% beta + 0.1*rnorm(n)
#' data = list(x = X, y = y)
#' fit = FSGL(data, index)

## This is our main function FSGL and it is a packaged function.
## Our main code is not in this function
FSGL <- function(data, index, maxit = 1000, thresh = 0.001, min.frac = 0.1,
                 nlam = 20, gamma = 0.8, standardize = TRUE, verbose = FALSE, step = 1,
                 reset = 10, alpha = 0.95, lambdas = NULL, UseUpperBound = FALSE)
{

  ## centering (and potentially scaling)
  out <- center_scale(data$x, standardize)
  data$x <- out$x
  X.transform <- out$X.transform

  ## linear model

  ## function FsglLinear is our main code for linear model calling C++ code after compilation
  Sol <- FsglLinear(data, index, thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, reset = reset, alpha = alpha, UseUpperBound=UseUpperBound)



  Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, X.transform = X.transform)

  class(Sol) = "FSGL"
  return(Sol)
}


## This function always centers X (and scales if standardize == TRUE)
## It also returns the centering/scaling numbers (for use in prediction with new data)
center_scale <- function(X, standardize){
  means <- apply(X,2,mean)
  X <- t(t(X) - means)

  X.transform <- list(X.means = means)

  if(standardize == TRUE){
    var <- apply(X,2,function(x)(sqrt(sum(x^2))))
    X <- t(t(X) / var)
    X.transform$X.scale <- var
  }
  else{
    X.transform$X.scale <- 1
  }

  return(list(x = X, X.transform = X.transform))
}
