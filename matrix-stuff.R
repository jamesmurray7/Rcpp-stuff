library(Rcpp)
library(RcppArmadillo)

setwd("~/Documents/Rcpp-stuff/")
sourceCpp("arma-matrix-experiments.cpp")

RexpdiagM <- function(x, M){
  crossprod(M, diag(exp(x)) %*% M)
}
RexpdiagM <- function(x, M){
  D <- diag(exp(x))
  crossprod(M, D %*% M)
}

x <- rnorm(10)
M <- matrix(rnorm(100), nrow = 10, ncol = 10)

bench <- microbenchmark::microbenchmark(
   R = RexpdiagM(x, M), R2 = RexpdiagMM(x, M), Cpp1 = expdiagM(x, M), Cpp2 = exdiagMM(x, M)
)

bench
plot(bench)


A <- B <- C <- D <- E <- matrix(rnorm(500*500),500,500)

bench <- microbenchmark::microbenchmark(
  R = A + B + C + D + E, Rreduce = Reduce('+', list(A,B,C,D,E)),
  Cpp = matrixSum(A,B,C,D,E)
)
bench
