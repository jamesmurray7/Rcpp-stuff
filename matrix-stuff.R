library(Rcpp)
library(RcppArmadillo)

setwd("~/Documents/Rcpp-stuff/")
sourceCpp("arma-matrix-experiments.cpp")

RexpdiagMM <- function(x, M){
  crossprod(M, diag(exp(x)) %*% M)
}

RexpdiagM <- function(x, M){
  D <- diag(exp(x))
  crossprod(M, D %*% M)
}

x <- rnorm(10)
M <- matrix(rnorm(100), nrow = 10, ncol = 10)

bench <- microbenchmark::microbenchmark(
   R = RexpdiagM(x, M), R2 = RexpdiagMM(x, M), 
   Cpp1 = expdiagM(x, M), Cpp2 = exdiagMM(x, M),
   times = 1e4
)

bench
plot(bench)


A <- B <- C <- D <- E <- matrix(rnorm(500*500),500,500)

bench <- microbenchmark::microbenchmark(
  R = A + B + C + D + E, Rreduce = Reduce('+', list(A,B,C,D,E)),
  Cpp = matrixSum(A,B,C,D,E)
)
bench
plot(bench)

# Crossprods
x <- rnorm(1000)
bench <- microbenchmark::microbenchmark(
  `R-tcrossprod` = {tcrossprod(x)},
  `R-outer` = {outer(x, x)},
  `C-tcrossprod` = {Ctcp(x)},
  times = 1e3
)
bench
plot(bench)
