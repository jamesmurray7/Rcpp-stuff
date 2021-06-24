#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//colvec bLL(vec& b, 
		   //const colvec& Y, const mat& X, const mat& Z, const mat& V, const mat& D,
		   //int mi, const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
		   //const rowvec& l0u, const mat& F, const vec& g, const vec& beta, const vec& eta,
		   //const vec& gr, const rowvec& rvFi){
    //colvec resid = Y - X * beta - Z * b;
    //Rcpp::Rcout << "resid: " << resid << std::endl;
    //double temp, obj;
    //if(Delta == 0){
		//double temp = 0;
	//}else{
		//double temp = log(l0i);
	//}
    //Rcpp::Rcout << "temp: " << temp << std::endl;
    //obj = -mi/2 * log(2 * datum::pi) - 1/2 * log(det(V)) - 1/2 * resid.t() * V.i() * resid +
		  //-5/2  * log(2 * datum::pi) - 1/2 * log(det(D)) - 1/2 * b.t() * D.i() * b;
    
	//return resid;
//}

//b.ll.nlm <- function(b.hat,			
                     //Y, X, Z, V, D, mi,			
                     //K, Delta, l0i, Fi, l0u, Fu, g, beta, eta){			
  //b <- c(b.hat[1], b.hat[2], b.hat[3], b.hat[4], b.hat[5], b.hat[6], b.hat[7], b.hat[8], b.hat[9], b.hat[10])			
  //gr <- rep(g, each = 2)			
  //if(Delta == 0) temp <- 0 else temp <- log(l0i)			
  //ll <- -5 * log(2*pi) - 1/2 * log(det(D)) - 1/2 * crossprod(b, solve(D) %*% b) - 			
    //mi/2 * log(2*pi) - 1/2 * log(det(V)) - 1/2*crossprod(Y - X %*% beta - Z %*% b, solve(V) %*% (Y - X %*% beta - Z%*%b)) + 			
    //temp + Delta * (K %*% eta + repVec(Fi) %*% (gr*b)) - l0u %*% (exp(K %*% eta) %x% exp(repCols(Fu) %*% (gr * b)))			
  //-ll			
//}

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
	double rtr = as_scalar(resid.t() * V.i() * resid);
	double C = -mi/2 * log(2 * M_PI);
	double ldV = as_scalar(log(det(V)));
	// Rcpp::Rcout << "resid: " << resid << std::endl;
	// Rcpp::Rcout << "RTVR: " << resid.t() * V.i() * resid << std::endl;
	return C - 1/2 * ldV -1/2 * rtr;
}

// [[Rcpp::export]]
double getRE(vec& b, mat& D){
	Rcpp::Rcout << "bDb: " << b.t() * D.i() * b;
	Rcpp::Rcout << "logDetD: " << log(det(D)) << std::endl;
	Rcpp::Rcout << "-5/2 * log(2 * pi)" << -5/2 * log(2 * M_PI) << std::endl;
	double C = -5/2 * log(2 * M_PI);
	double lDD = as_scalar(log(det(D)));
	double bDb = as_scalar(b.t() * D.i() * b);
	return C - 1/2 * lDD - 1/2 * bDb;
}

// [[Rcpp::export]]
mat repMat(const mat& M, int K){
	return repmat(M, 1, K);
}

