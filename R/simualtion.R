# # # ##simulation
# # #
# # library(microbenchmark)
# # library(SGL)
# # library(FSGL)
# # # source('./R/betterPathCalc.R')
# # # source('./R/FsglLinear.R')
# # # source('./R/FSGL.R')
# # # Rcpp::sourceCpp('src/fsgl.cpp')
# # # options(scipen = 200)
# # # options(digits=5)
# n = 100; p =20; size.groups = 5
# index <- ceiling(1:p / size.groups)
# set.seed(545564)
# X = matrix(rnorm(n * p), ncol = p, nrow = n)
# beta = (-2:2)
# y = X[,1:5] %*% beta + 0.1*rnorm(n)
# data = list(x = X, y = y)
# # #
# # # # microbenchmark::microbenchmark(
# # # #   FSGL(data, index, alpha = 0.85, UseUpperBound = FALSE),
# # # #   SGL(data, index, alpha = 0.85),
# # # #   times = 1
# # # # )
# # # t0 = Sys.time()
# fit1_simu = cvFSGL(data, index, lambdas=c(0.21,0.11),UseUpperBound = FALSE)
# # # t1 = Sys.time()
# # # fit2_simu = SGL(data, index, type = "linear")
# # # t2 = Sys.time()
# # #
# # # t1-t0
# # # t2-t1
# # # # #
# # #
# # # # cat("the same number of beta between FSGL and SGL is", sum(fit1$beta==fit2$beta),
# # # #     "the total number of beta is", nrow(fit1$beta)*ncol(fit1$beta),
# # # #     ", so most beta of SGL and FSGL are the same ")
# # #
# # #
# # # ## real data
# # #
# # abalone <- read.csv("./data/abalone.csv", header=FALSE)
# # # # # # ### use v2 to v8 as variable
# # # # # # ### each pair (a, b) to generate 5-dimensional vector (a,b,a*b,a^2,b^2)
# # # abalone <- with(abalone,{
# # #   adalone$V1 <- factor(abalone$V1)
# # # })
# # size.groups = choose(7,2)
# # p = size.groups * 5
# # index = ceiling(1:p / size.groups)
# # X = matrix(0, nrow(abalone), p)
# # count = 0
# # ### create group index and beta
# # for(i in 1:7){
# #   for(j in 1:i){
# #     if(i != j){
# #       X[,count+1] = abalone[,i+1]
# #       X[,count+2] = abalone[,j+1]
# #       X[,count+3] = abalone[,i+1] * abalone[,j+1]
# #       X[,count+4] = abalone[,i+1]^2
# #       X[,count+5] = abalone[,j+1]^2
# #       count = count+5
# #     }
# #   }
# # }
# # Y = abalone[,9]
# # data = list(x = X, y = Y)
# # t01 = Sys.time()
# # fit1_real = FSGL(data, index, alpha=0.95)
# # t11 = Sys.time()
# # fit2_real = SGL(data, index, alpha=0.95)
# # t21 = Sys.time()
# #
# # t11-t01
# # t21-t11
# # # # # cat("the same number of beta between FSGL and SGL is", sum(fit1_real$beta==fit2_real$beta),
# # # # #     "the total number of beta is", nrow(fit1_real$beta)*ncol(fit1_real$beta),
# # # # #     ", so most beta of SGL and FSGL are the same \n")
# # # # # #
# # # # delta1 = t1-t0
# # # # delta2 = t2-t1
# # # # cat("the running time of FSGL is", delta1,", running time of SGL is", delta2,"\n")
# # # # mse = sum((X%*%as.matrix(apply(fit2_real$beta,1,mean)) - Y)^2)
# # # # cat("the mse is", mse)
# # # # # # microbenchmark::microbenchmark(
# # # # # #   FSGL(data, index, alpha = 0.6),
# # # # # #   SGL(data, index, alpha = 0.6),
# # # # # #   times = 100
# # # # # # )
# # # #
# # # #
# # bodyfat <- read.table("./data/bodyfat.csv", header=FALSE, sep=',')
# # ### use v2 to v8 as variable
# # ### each pair (a, b) to generate 5-dimensional vector (a,b,a*b,a^2,b^2)
# # size.groups = choose(14,2)
# # p = size.groups * 5
# # index = ceiling(1:p / size.groups)
# # X = matrix(0, nrow(bodyfat), p)
# # count = 0
# # ### create group index and beta
# # for(i in 1:14){
# #   for(j in 1:i){
# #     if(i != j){
# #       X[,count+1] = bodyfat[,i+1]
# #       X[,count+2] = bodyfat[,j+1]
# #       X[,count+3] = bodyfat[,i+1] * bodyfat[,j+1]
# #       X[,count+4] = bodyfat[,i+1]^2
# #       X[,count+5] = bodyfat[,j+1]^2
# #       count = count+5
# #     }
# #   }
# # }
# # Y = bodyfat[,15]
# # data = list(x = X, y = Y)
# # t02 = Sys.time()
# # fit1_real = FSGL(data, index, alpha=0.95)
# # t12 = Sys.time()
# # fit2_real = SGL(data, index, alpha=0.95)
# # t22 = Sys.time()
# # t12-t02
# # t22-t12
# # # # cat("the same number of beta between FSGL and SGL is", sum(fit1_real$beta==fit2_real$beta),
# # # #     "the total number of beta is", nrow(fit1_real$beta)*ncol(fit1_real$beta),
# # # #     ", so most beta of SGL and FSGL are the same \n")
# # # #
# # # # cat("the running time of FSGL is", t1-t0,", running time of SGL is", t2-t1,"\n")
# # #
# # # #
# # white<- read.csv("./data/white.csv", header=TRUE)
# # ### use v2 to v8 as variable
# # ### each pair (a, b) to generate 5-dimensional vector (a,b,a*b,a^2,b^2)
# # size.groups = choose(11,2)
# # p = size.groups * 5
# # index = ceiling(1:p / size.groups)
# # X = matrix(0, nrow(white), p)
# # count = 0
# # ### create group index and beta
# # for(i in 1:11){
# #   for(j in 1:i){
# #     if(i != j){
# #       X[,count+1] = white[,i]
# #       X[,count+2] = white[,j]
# #       X[,count+3] = white[,i] * white[,j]
# #       X[,count+4] = white[,i]^2
# #       X[,count+5] = white[,j]^2
# #       count = count+5
# #     }
# #   }
# # }
# # Y = white[,12]
# # data = list(x = X, y = Y)
# # t03 = Sys.time()
# # fit1 = FSGL(data, index,alpha=0.5)
# # t13 = Sys.time()
# # fit2 = SGL(data, index,alpha=0.5)
# # t23 = Sys.time()
# # t13-t03
# # t23-t13
# # # # # delta1=t1-t0
# # # # # delta2=t2-t1
# # # # # # cat("the same number of beta between FSGL and SGL is", sum(fit1$beta==fit2$beta),
# # # # #     "the total number of beta is", nrow(fit1$beta)*ncol(fit1$beta),
# # # # #     ", so most beta of SGL and FSGL are the same \n")
# # # #
# # # # # cat("the running time of FSGL is",delta1,", running time of SGL is", delta2,"\n")
# # # # t1-t0
# # # # t2-t1
# # # # mse = sum((X%*%as.matrix(apply(fit2$beta,1,mean)) - Y)^2)
# # # # cat("the mse is", mse)
# # # #
# # # #
# # # # ## red
# # # #
# # # # red<- read.csv("./data/red.csv", header=TRUE)
# # # # # ### use v2 to v8 as variable
# # # # # ### each pair (a, b) to generate 5-dimensional vector (a,b,a*b,a^2,b^2)
# # # # size.groups = choose(11,2)
# # # # p = size.groups * 5
# # # # index = ceiling(1:p / size.groups)
# # # # X = matrix(0, nrow(red), p)
# # # # count = 0
# # # # ### create group index and beta
# # # # for(i in 1:11){
# # # #   for(j in 1:i){
# # # #     if(i != j){
# # # #       X[,count+1] = red[,i]
# # # #       X[,count+2] = red[,j]
# # # #       X[,count+3] = red[,i] * red[,j]
# # # #       X[,count+4] = red[,i]^2
# # # #       X[,count+5] = red[,j]^2
# # # #       count = count+5
# # # #     }
# # # #   }
# # # # }
# # # # Y = red[,12]
# # # # data = list(x = X, y = Y)
# # # # t0 = Sys.time()
# # # # fit1 = FSGL(data, index,alpha=0.8)
# # # # t1 = Sys.time()
# # # # fit2 = SGL(data, index,alpha=0.8)
# # # # t2 = Sys.time()
# # # # # # # delta1=t1-t0
# # # # # # # delta2=t2-t1
# # # # # #
# # # # mse = sum((X%*%as.matrix(apply(fit2$beta,1,mean)) - Y)^2)
# # # # cat("the mse is", mse)
# # # # t1-t0
# # # # t2-t1
# # # # # cat("the running time of FSGL is",deltat1,", running time of SGL is", t2-t1,"\n")
# # #
