cross.terms <- matrix(0,5,5)
for(L in 1:4){
  for(M in 2:5){
    if(L!=M){
      for(k in 1:gh.nodes){
        xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
        cross.terms[M,L] <- cross.terms[M,L] + 
          w[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * Fu[[i]] %*% bb[M, ]) + 
          gamma[M] * v[k] * w[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * tau2.surv * tau.tilde[M, ]) + 
          2 * gamma[L] * v[k] * w[k] * crossprod(xi * tau.surv * xi, Fu[[i]] %*% bb[M, ]) + 
          2 * gamma[L] * gamma[M] * w[k] * v[k]^2 * crossprod(xi * tau.surv * xi * tau2.surv, tau.tilde[M, ]) + 
          gamma[L] * gamma[M] * v[k] * w[k] * crossprod(xi * xi * tau2.surv, tau.tilde[M, ])
      }
    } 
  }
}
library(Rcpp)
library(RcppArmadillo)
sourceCpp('~/Documents/Rcpp-stuff/gamma-calc.cpp')
# function (gamma, tautilde, tausurv, tau2surv, musurv, w, v, Fu, 
#           haz, bb, L, gh) 
gammaCalc(gamma, tau.tilde, tau.surv, tau2.surv, mu.surv, w, v, Fu[[i]], l0u[[i]], bb, L = 5, gh = 3)


bench <- microbenchmark::microbenchmark(
  R = {Igamma.store <- matrix(NA, nr = gh.nodes, nc = 5)
  # Cross-terms followed by second derivatives //
  gammacalcstore[[i]] <- gammaCalc(gamma, tau.tilde, tau.surv, tau2.surv, mu.surv, w, v, Fu[[i]], l0u[[i]], bb, 5, gh.nodes)
  cross.terms <- matrix(0,5,5)
  for(L in 1:4){
    for(M in 2:5){
      if(L!=M){
        for(k in 1:gh.nodes){
          xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
          cross.terms[M,L] <- cross.terms[M,L] + 
            w[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * Fu[[i]] %*% bb[M, ]) + 
            gamma[M] * v[k] * w[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * tau2.surv * tau.tilde[M, ]) + 
            2 * gamma[L] * v[k] * w[k] * crossprod(xi * tau.surv * xi, Fu[[i]] %*% bb[M, ]) + 
            2 * gamma[L] * gamma[M] * w[k] * v[k]^2 * crossprod(xi * tau.surv * xi * tau2.surv, tau.tilde[M, ]) + 
            gamma[L] * gamma[M] * v[k] * w[k] * crossprod(xi * xi * tau2.surv, tau.tilde[M, ])
        }
      } 
    }
  }
  for(k in 1:gh.nodes){
    xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
    for(L in 1:5){
      Igamma.store[k, L] <- w[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * Fu[[i]] %*% bb[L, ]) + 
        gamma[L] * w[k] * v[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * tau2.surv * tau.tilde[L, ]) + 
        v[k] * w[k] * crossprod(xi * tau.surv, xi) + 
        2 * gamma[L] * v[k] * w[k] * crossprod(xi * tau.surv * xi, Fu[[i]] %*% bb[L, ]) + 
        2 * gamma[L]^2 * v[k]^2 * w[k] * crossprod(xi * tau.surv * xi * tau2.surv, tau.tilde[L, ]) + 
        gamma[L]^2 * v[k] * w[k] * crossprod(xi * xi * tau2.surv, tau.tilde[L, ])
    }
  }
  
  Igammai <- colSums(Igamma.store)
  Igammai <- diag(an(Igammai))
  Igammai[lower.tri(Igammai)] <- cross.terms[lower.tri(cross.terms)]
  Igammai[upper.tri(Igammai)] <- t(cross.terms)[upper.tri(t(cross.terms))]},
  Cpp = gammaCalc(gamma, tau.tilde, tau.surv, tau2.surv, mu.surv, w, v, Fu[[i]], l0u[[i]], bb, L = 5, gh = 3)
)

Rtimes <- bench[bench$expr=="R", 2]
Rtimes <- Rtimes[which(Rtimes < quantile(Rtimes, .95))]
Cpptimes <- bench[bench$expr!="R", 2]
Cpptimes <- Cpptimes[which(Cpptimes < quantile(Cpptimes, .95))]

cbind(R = Rtimes, Cpp = Cpptimes) %>% 
  as.data.frame %>% 
  gather() %>% 
  ggplot(aes(y = value, x = key)) + 
  geom_boxplot() + 
  scale_y_log10()
  theme_minimal()
