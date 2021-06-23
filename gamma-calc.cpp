#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//# I(\gamma) ----
      //Igamma.store <- matrix(NA, nr = gh.nodes, nc = 5)
      //# Cross-terms followed by second derivatives //
      //cross.terms <- matrix(0,5,5)
      //for(L in 1:4){
        //for(M in 2:5){
          //if(L!=M){
            //for(k in 1:gh.nodes){
              //xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
              //cross.terms[M,L] <- cross.terms[M,L] + 
                //w[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * Fu[[i]] %*% bb[M, ]) + 
                //gamma[M] * v[k] * w[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * tau2.surv * tau.tilde[M, ]) + 
                //2 * gamma[L] * v[k] * w[k] * crossprod(xi * tau.surv * xi, Fu[[i]] %*% bb[M, ]) + 
                //2 * gamma[L] * gamma[M] * w[k] * v[k]^2 * crossprod(xi * tau.surv * xi * tau2.surv, tau.tilde[M, ]) + 
                //gamma[L] * gamma[M] * v[k] * w[k] * crossprod(xi * xi * tau2.surv, tau.tilde[M, ])
            //}
          //} 
        //}
      //}
      //for(k in 1:gh.nodes){
        //xi <- l0u[[i]] * (mu.surv * exp(v[k] * tau.surv))
        //for(L in 1:5){
          //Igamma.store[k, L] <- w[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * Fu[[i]] %*% bb[L, ]) + 
            //gamma[L] * w[k] * v[k] * t(bb[L, ]) %*% t(Fu[[i]]) %*% (xi * tau2.surv * tau.tilde[L, ]) + 
            //v[k] * w[k] * crossprod(xi * tau.surv, xi) + 
            //2 * gamma[L] * v[k] * w[k] * crossprod(xi * tau.surv * xi, Fu[[i]] %*% bb[L, ]) + 
            //2 * gamma[L]^2 * v[k]^2 * w[k] * crossprod(xi * tau.surv * xi * tau2.surv, tau.tilde[L, ]) + 
            //gamma[L]^2 * v[k] * w[k] * crossprod(xi * xi * tau2.surv, tau.tilde[L, ])
        //}
      //}

// So we need 
// Vectors: gamma, tau.surv, tau2.surv, mu.surv, w, v, Fu, l0u.
// Matrix: bb (matrix of bi's); tau.tilde (matrix STORE of diag vecs)

// i = ROWS (M)
// j = COLS (L)


// Quick test re: col/rowvecs
// [[Rcpp::export]]
rowvec qtest(mat tautilde, int i){
  arma::rowvec a = tautilde.row(i);
  return a;
}

// [[Rcpp::export]]
mat gammaCalc(const vec& gamma,const mat& tautilde, const vec& tausurv, const vec& tau2surv, const colvec& musurv, 
              const vec& w, const vec& v, const vec& Fu, const vec& haz, const mat& bb, int L, int gh){
	mat M = zeros<mat>(L, L);
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			const rowvec bL = bb.row(j);
			const rowvec bM = bb.row(i);
			if(i!=j){
				for(int k = 0; k < gh; k++){
					const colvec xi = haz % (musurv * exp(v(k) * tausurv));
					const colvec temp1 = xi % tausurv % xi;
					const colvec temp2 = temp1 % tau2surv;
					const colvec temp3 = xi % xi % tau2surv;
					M(i,j) += as_scalar(w(k) * bL.t() * Fu.t() * (xi % Fu * bM) + 
					          gamma[i] * w[k] * v[k] * bL.t() * Fu.t() * (xi % tau2surv % tautilde.row(i)) +
					          2 * gamma[j] * w[k] * v[k] * temp1.t() * (Fu * bM) + 
					          2 * gamma[i] * gamma[j] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i) + 
					          gamma[i] * gamma[j] * v[k] * w[k] * temp3.t() * tautilde.row(i));
				}
			}
		}
	}
	return M;
}

// Testing JUST the creation of xi with this ---
// [[Rcpp::export]]
colvec vecsum(const rowvec& tausurv, const colvec& musurv, 
              vec& v, const rowvec& haz, int gh){
  rowvec xi = zeros<rowvec>(30);
  Rcpp::Rcout << "haz: " << haz << std::endl;
  Rcpp::Rcout << "musurv: " << musurv << std::endl;
	for(int k = 0; k < gh; k++){
		Rcpp::Rcout << "v:" << v(k) << std::endl;;
		Rcpp::Rcout << "v * tau.surv:" << v(k) * tausurv << std::endl;
		Rcpp::Rcout << "musurv * tau.surv:" << musurv.t() % tausurv << std::endl;
	  xi += haz % (musurv.t() % exp(v(k) * tausurv));
	}
	return xi.t();
}

