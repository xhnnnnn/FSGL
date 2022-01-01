FsglLinear <-
  function(data, index, thresh = 0.0001, nlam = 20, lambdas = NA, inner.iter = 100, outer.iter = 100, outer.thresh = 0.0001, gamma = 0.8, step = 1, reset = 10, alpha = 0.95, min.frac = 0.05, verbose = FALSE, UseUpperBound=FALSE){
    ## If not given lambdas, we use a method to generate lambda
    if(exists("lambdas")){
      lambdas <- NA;
    }
    if(is.na(lambdas)){
      lambdas <- betterPathCalc(data = data, index = index, alpha=alpha, min.frac = min.frac, nlam = nlam, type = "linear")
    }

    X <- data$x
    y <- data$y
    n <- nrow(X)
    p <- ncol(X)

    ## Setting up group lasso stuff

    ord <- order(index)
    index <- index[ord]
    X <- X[,ord]
    unOrd <- match(1:length(ord),ord)

    ## setting up group index and group range index

    groups <- unique(index)
    num.groups <- length(groups)
    range.group.ind <- rep(0,(num.groups+1))
    for(i in 1:num.groups){
      range.group.ind[i] <- min(which(index == groups[i])) - 1
    }
    range.group.ind[num.groups+1] <- ncol(X)

    group.length <- diff(range.group.ind)

    ## Calling C++ code to run the algorithms using different lambda

    #alpha <- sqrt(2*log(p))/(1+sqrt(2*log(num.groups)/min(group.length)) + sqrt(2*log(p)))
    nlam = length(lambdas)
    # beta.is.zero <-array(0, c(ncol(X),nlam))
    beta <- array(0, c(ncol(X),nlam))

    # eta <- array(0, c(n, nlam))7

    junk <- FastSGL(X = X, y = y, numGroup = num.groups, rangeGroupInd =  range.group.ind, groupLen = group.length, alpha=alpha, lambdas = lambdas, innerIter = inner.iter, outerIter = outer.iter, thresh = as.double(thresh), outerThresh = as.double(outer.thresh), gamma = as.double(gamma), step = as.double(step), reset = as.integer(reset), UseUPB = as.numeric(UseUpperBound))
    beta <- junk$beta
    # eta <- junk$eta
    # beta.is.zero <- junk$betaIsZero

    if(verbose==TRUE){
      for(k in lambdas){
        write(paste("***Lambda", k, "***"),"")
      }

    }


    return(list(beta = beta[unOrd,], lambdas = lambdas))
  }

