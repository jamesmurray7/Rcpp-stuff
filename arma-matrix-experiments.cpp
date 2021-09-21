#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Doing random experiments
// First two are two competing ways of calculating 
// D <- diag(exp(v))
// return(crossprod(M, D %*% M)) where v vector and M matrix

// [[Rcpp::export]]
mat expdiagM(const vec& x, const mat& M){
	colvec diagonal(x.size());
	for(int i = 0; i < x.size(); i++){
		diagonal(i) = exp(x[i]);
	}
	mat D = diagmat(diagonal);
	return M.t() * D * M;
}

// And another way of doing the same thing, this should be faster
// As we're telling Rcpp what diagvec is 'upfront'.

// [[Rcpp::export]]
mat exdiagMM(const vec& x, const mat& M){
	return M.t() * diagmat(exp(x)) * M;
}

// Simply adding together matrices

// [[Rcpp::export]]
mat matrixSum(mat& A, mat& B, mat& C, mat& D, mat& E){
	return A + B + C + D + E;
}

// Generate empty matrix of size L

// [[Rcpp::export]]
mat genZeros(int& L){
	mat M = zeros<mat>(L,L);
	return M;
}

// Generate empty matrix and then loop through elements, adding i+j
// [[Rcpp::export]]
mat loopMat(int& L){
	mat M = zeros<mat>(L,L);
	for(int i = 1; i <= L; i++){
		for(int j = 1; j <= L; j++){
		//	Rcpp::Rcout << "i:" << i;
		//	Rcpp::Rcout << "j:" << j << std::endl;
			M(i-1,j-1) += i*j;
		}
	}
	return M;
}

// [[Rcpp::export]]
mat Ctcp(vec& x){
	return x * x.t();
}

