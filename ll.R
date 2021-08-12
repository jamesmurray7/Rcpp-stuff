library(Rcpp)
library(RcppArmadillo)

sourceCpp("../../gammaCalc.cpp")
sourceCpp("../../../Rcpp-stuff/ll.cpp")

# longitudinal
longit <- function(b, Y, X, Z, beta, V, mi){
  resid <- Y - X %*% beta - Z %*% b
  return(-mi/2 * log(2 * pi) - 1/2 * log(det(V)) - 1/2 * crossprod(resid, solve(V) %*% resid))
}

longit(bb[1,], Y[[1]], X[[1]], Z[[1]], beta, V[[1]], sum(mi[[1]]))
getLongit(bb[1,], Y[[1]], X[[1]], Z[[1]], beta, V[[1]], sum(mi[[1]]))

# REs
res <- function(b, D){
  -5 * log(2 * pi) - 1/2 * log(det(D)) - 1/2 * crossprod(b, solve(D) %*% b)
}

res(bb[1,], D)
getRE(bb[1,], D)  

#survival
surv <- function(b, K, Delta, l0i, Fi, l0u, Fu, g, beta, eta){
  gr = rep(gamma, each = 2)
  temp <- 0
  if(Delta == 1) temp <- log(l0i)
  temp + Delta * (K %*% eta + repVec(Fi) %*% (gr*b)) - l0u %*% (exp(K %*% eta) %x% exp(repCols(Fu) %*% (gr * b)))
}

surv(bb[1,], K[[1]], Di[1], l0i[1], Fi[1,], l0u[[1]], Fu[[1]], gamma, beta, eta)
getSurvival(bb[1,], K[[1]], Di[1], l0i[1], Fi[1,], l0u[[1]], Fu[[1]], gamma, beta, eta, rep(gamma, each = 2), repVec(Fi[1,]))

# log likelihood
b.ll.nlm(bb[1,], Y[[1]], X[[1]], Z[[1]], V[[1]], D, sum(mi[[1]]), K[[1]], Di[1],
         l0i[1], Fi[1,], l0u[[1]], Fu[[1]], gamma, beta, eta)
bllC(bb[1,], Y[[1]], X[[1]], Z[[1]], V[[1]], D, sum(mi[[1]]), K[[1]], Di[1], 
     l0i[1], Fi[1,], l0u[[1]], Fu[[1]], gamma, beta, eta, rep(gamma, each = 2), repVec(Fi[1,]))

# gradient 
nlminb.gradient(bb[1,], Y[[1]], X[[1]], Z[[1]], V[[1]], D, sum(mi[[1]]), K[[1]], Di[1],
                l0i[1], Fi[1,], l0u[[1]], Fu[[1]], gamma, beta, eta)
gradC(bb[1,], Y[[1]], X[[1]], Z[[1]], V[[1]], D, sum(mi[[1]]), K[[1]], Di[1], 
     l0i[1], Fi[1,], l0u[[1]], Fu[[1]], gamma, beta, eta, rep(gamma, each = 2), repVec(Fi[1,]))

microbenchmark::microbenchmark(
  R = {
    llR <- c()
    for(i in 1:100){
      llR[i] <- b.ll.nlm(bb[i,], Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
                        l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta)
    }
  },
  C = {
    llC <- c()
    for(i in 1:100){
      llC[i] <- bllC(bb[i,], Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
                        l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta, rep(gamma, each = 2), repVec(Fi[1,]))
    }
  }, times = 100
) -> bench


nlm(b.ll.nlm, bb[i,], Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
    l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta)
nlm(bllC, bb[i,], Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
    l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta, rep(gamma, each = 2), repVec(Fi[1,]))

bench2 <- microbenchmark::microbenchmark(
  R_nlm = {
    bbR <- matrix(NA, 100, 10)
    for(i in 1:100){
      bbR[i, ] <- nlm(b.ll.nlm, bb[i,], Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
                      l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta)$estimate
    }
  },
  C_nlm = {
    bbC <- matrix(NA, 100, 10)
    for(i in 1:100){
      bbC[i, ] <- nlm(bllC, bb[i,], Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
                      l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta, rep(gamma, each = 2), repVec(Fi[1,]))$estimate
    }
  }, times = 10
)

bench3 <- microbenchmark::microbenchmark(
  R_ucminf = {
    bbR <- matrix(NA, 100, 10)
    for(i in 1:100){
      bbR[i,] <- ucminf::ucminf(bb[i, ], b.ll.nlm, nlminb.gradient, 
                                Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
                                l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta,
                                hessian = 0)$par
    }
  },
  R_ucminf_steptol = {
    bbRs <- matrix(NA, 100, 10)
    for(i in 1:100){
      bbRs[i,] <- ucminf::ucminf(bb[i, ], b.ll.nlm, nlminb.gradient, 
                                Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
                                l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta,
                                hessian = 0, control = list(grtol = 1e-3, xtol = 1e-3))$par
    }
  },
  C_ucminf = {
    bbC <- matrix(NA, 100, 10)
    for(i in 1:100){
      bbC[i, ] <- ucminf::ucminf(bb[i,], bllC, gradC,
                                 Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
                                 l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta, rep(gamma, each = 2), repVec(Fi[1,]),
                                 hessian= 0)$par
    }
  },
  C_ucminf_steptol = {
    bbCs <- matrix(NA, 100, 10)
    for(i in 1:100){
      bbCs[i, ] <- ucminf::ucminf(bb[i,], bllC, gradC,
                                 Y[[i]], X[[i]], Z[[i]], V[[i]], D, sum(mi[[i]]), K[[i]], Di[i],
                                 l0i[i], Fi[i,], l0u[[i]], Fu[[i]], gamma, beta, eta, rep(gamma, each = 2), repVec(Fi[1,]),
                                 hessian = 0, control = list(grtol = 1e-3, xtol = 1e-3))$par
    }
  },
  times = 50
)

bench3
plot(bench3)
mean(abs(bbC-bbR))
