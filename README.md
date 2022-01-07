# FSGL

R package of Fast Sparse Group Lasso (beta version)

Installation
------------

`FSGL` can be installed via Github using `devtools`

    # install.packages("devtools")
    library(devtools)
    devtools::install_github("xhnnnn/FSGL")

You'll need a working C++11 compiler, which can obtained by installing
Xcode on MacOS, and RTools on Windows.

Example
-------

    library(FSGL)
    
    n = 100; p = 20; size.groups = 10
    index <- ceiling(1:p / size.groups)
    X = matrix(rnorm(n * p), ncol = p, nrow = n)
    beta = (-2:2)
    y = X[,1:5] %*% beta + 0.1*rnorm(n)
    data = list(x = X, y = y)
    fit = FSGL(data, index)



References
----------

Ida Y, Fujiwara Y, Kashima H. Fast Sparse Group Lasso[J]. Advances in Neural Information Processing Systems, 2019, 32: 1702-1710.