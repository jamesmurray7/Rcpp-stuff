#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ============================
// INDIVIDUAL TESTING FUNCTIONS 
// ============================

// [[Rcpp::export]]
colvec getResid(vec& b, const colvec& Y, const mat& X, const mat& Z, const vec& beta){
	colvec resid = Y - X * beta - Z * b;
	return resid;
}

// [[Rcpp::export]]
double getLongit(vec& b, const colvec& Y, const mat& X, const mat& Z, const vec& beta, const mat& V, int mi){
	colvec resid = Y - X * beta - Z * b;
	return as_scalar(-mi/2.0 * log(2.0 * M_PI) - 0.5 * log(det(V)) -0.5 * resid.t() * V.i() * resid);
}

// [[Rcpp::export]]
double getRE(vec& b, mat& D){
	return as_scalar(-5.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b);
}

// [[Rcpp::export]]
mat repMat(const mat& M, int K){
	return repmat(M, 1, K);
}

// [[Rcpp::export]]
double getSurvival(vec& b, const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
                   const rowvec& l0u, const mat& Fu, const vec& g, const vec& beta, const vec& eta,
				   const vec& gr, const rowvec& rvFi){
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return as_scalar(temp + Delta * (K * eta + rvFi * (gr % b)) - l0u * (kron(exp(K * eta), exp(repmat(Fu, 1, 5) * (gr % b)))));
}

// ============================
// ALL TOGETHER
// ============================

// [[Rcpp::export]]
double bllC(vec& b, const colvec& Y, const mat& X, const mat& Z, const mat& V, const mat& D,
		   int mi, const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
		   const rowvec& l0u, const mat& Fu, const vec& g, const vec& beta, const vec& eta,
		   const vec& gr, const rowvec& rvFi){
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	colvec resid = Y - X * beta - Z * b;
	return -1.0 * as_scalar(-mi/2.0 * log(2.0 * M_PI) - 0.5 * log(det(V)) -0.5 * resid.t() * V.i() * resid +
	                 -5.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
	                 temp + Delta * (K * eta + rvFi * (gr % b)) - l0u * (kron(exp(K * eta), exp(repmat(Fu, 1, 5) * (gr % b)))));
}

// Gradient function
// [[Rcpp::export]]
colvec gradC(vec& b, const colvec& Y, const mat& X, const mat& Z, const mat& V, const mat& D,
		     int mi, const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
		     const rowvec& l0u, const mat& Fu, const vec& g, const vec& beta, const vec& eta,
		     const vec& gr, const rowvec& rvFi){
	colvec resid = Y - X * beta - Z * b;
	return -1.0 * (Z.t() * V.i() * resid - D.i() * b + Delta * rvFi.t() % gr + 
	                 - repmat(Fu, 1, 5).t() * (l0u.t() % (kron(exp(K * eta), exp(repmat(Fu, 1, 5) * (gr % b))))) % gr);

}

// function(b.hat,
//                             Y, X, Z, V, D, mi,
//                             K, Delta, l0i, Fi, l0u, Fu, g, beta, eta){
//   b <- c(b.hat[1], b.hat[2], b.hat[3], b.hat[4], b.hat[5], b.hat[6], b.hat[7], b.hat[8], b.hat[9], b.hat[10])
//   gr <- rep(g, each = 2)
//   if(nrow(Fu) == 1){
//     db.ll <- crossprod(Z, solve(V) %*% (Y - X %*% beta - Z %*% b))  - solve(D) %*% b + 
//       Delta * repVec(Fi) * gr - repCols(Fu) %*% (l0u * (exp(K %*% eta) %x% exp(repCols(Fu) %*% (gr * b)))) * gr
//   }else{
//     db.ll <- crossprod(Z, solve(V) %*% (Y - X %*% beta - Z %*% b))  - solve(D) %*% b + 
//       Delta * repVec(Fi) * gr - crossprod(repCols(Fu), l0u * (exp(K %*% eta) %x% exp(repCols(Fu) %*% (gr * b)))) * gr
//   }
//   -db.ll
// }


