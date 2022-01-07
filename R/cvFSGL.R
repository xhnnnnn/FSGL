#' @title  Fit and Cross-Validate a Linear Model with a Combination of Lasso and Group Lasso Regularization
#'
#' @param data data should be a list of $x$ and $y$, x is a data matrix (n x p) and y is a vector
#' @param index A p-vector indicating group membership of each covariate
#' @param nfold Number of folds of the cross-validation loop
#' @param nlam Number of lambda to use in the regularization path
#' @param min.frac The minimum value of the penalty parameter, as a fraction of the maximum value
#' @param alpha The mixing parameter. \code{alpha = 1} is the lasso penalty.
#' @param lambdas  A user inputted sequence of lambda values for fitting. We recommend leaving this NULL and letting FSGL self-select values
#' @param thresh Convergence threshold for change in beta
#' @param maxit Maximum number of iterations to convergence
#' @param gamma Fitting parameter used for tuning backtracking (between 0 and 1)
#' @param verbose Logical flag for whether or not step number will be output
#' @param step Fitting parameter used for inital backtracking step size (between 0 and 1)
#' @param reset Fitting parameter used for taking advantage of local strong convexity in nesterov momentum (number of iterations before momentum term is reset)
#' @param foldid An optional user-pecified vector indicating the cross-validation fold in which each observation should be included. Values in this vector should range from 1 to nfold. If left unspecified, SGL will randomly assign observations to folds
#' @param UseUpperBound A logical flag for using upper bound
#'
#' @return
#' An object with S3 class "cv.FSGL"
#' \describe{
#'      \item{lldiff}{An \code{nlam} vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path)}
#'      \item{llSD}{An \code{nlame} vector of approximate standard deviations of \code{lldiff}}
#'      \item{lambdas}{The actual list of \code{lambda }values used in the regularization path.}
#'      \item{fit}{A model fit object created by a call to \code{FSGL} on the entire dataset}
#'      \item{foldid}{A vector indicating the cross-validation folds that each observation is assigned to}
#'      \item{prevals}{A matrix of prevalidated predictions for each observation, for each lambda-value}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' n = 50; p = 100; size.groups = 10
#' index <- ceiling(1:p / size.groups)
#' X = matrix(rnorm(n * p), ncol = p, nrow = n)
#' beta = (-2:2)
#' y = X[,1:5] %*% beta + 0.1*rnorm(n)
#' data = list(x = X, y = y)
#' cvFit = cvFSGL(data, index)
cvFSGL<-
  function(data, index, nfold = 10, nlam=20, min.frac = 0.05, alpha = 0.95,lambdas = NULL, thresh = 0.0001, maxit = 10000, gamma = 0.8, verbose = TRUE, step = 1, reset = 10, foldid = NA, UseUpperBound = FALSE){

    X <- data$x
    y <- data$y
    n <- nrow(X)
    p <- ncol(X)

    ## Setting up group lasso stuff ##

    ord <- order(index)
    index <- index[ord]
    X <- X[,ord]
    unOrd <- match(1:length(ord),ord)

    ## Coming up with other C++ info ##

    groups <- unique(index)
    num.groups <- length(groups)
    range.group.ind <- rep(0,(num.groups+1))
    for(i in 1:num.groups){
      range.group.ind[i] <- min(which(index == groups[i])) - 1
    }
    range.group.ind[num.groups+1] <- ncol(X)

    group.length <- diff(range.group.ind)
    beta.naught <- rep(0,ncol(X))
    beta <- beta.naught

    ## Done with group stuff ##

    ## finding the path

    MainSol <- FSGL(data, index, thresh = thresh, maxit = maxit, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, step = step, reset = reset, alpha = alpha, UseUpperBound = UseUpperBound)

    lambdas <- MainSol$lambdas
    nlam <- length(lambdas)
    lldiff <- rep(0, nlam)
    lldiffFold <- matrix(0, nrow = nlam, ncol = nfold)

    prevals <- matrix(0, nrow = nrow(data$x), ncol = nlam)
    if(sum(is.na(foldid))!=0){
      foldid = ceiling(1:n / nfold)
    }
    for(i in 1:nfold){
      ind.out <- which(foldid == i)
      ind.in <- which(foldid != i)
      new.data <- list(x = data$x[ind.in,], y = data$y[ind.in])

      new.sol <- FSGL(new.data, index, thresh = thresh, maxit = maxit, lambdas = lambdas, min.frac = min.frac, nlam = nlam, gamma = gamma, step = step, reset = reset, alpha = alpha, UseUpperBound=UseUpperBound)

      for(k in 1:nlam){
        etas <- X[ind.out,] %*% new.sol$beta[ord,k] ## Have to reorder betas according to ordering of index!
        lldiffFold[k,i] <- sum((y[ind.out] - etas)^2) / 2
        prevals[ind.out,k] <- etas
      }
      if(verbose == TRUE){
        write(paste("*** NFOLD ", i, "***"),"")
      }
    }
    lldiff = rowSums(lldiffFold)
    lldiffSD <- apply(lldiffFold,1,sd) * sqrt(nfold)
    obj <- list(lambdas = lambdas, lldiff = lldiff,llSD = lldiffSD, fit = MainSol, prevals = prevals)
    class(obj)="cv.SGL"
    return(obj)
  }
