#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// S is this -->
// lapply(Sigmai.store, function(x) lapply(split(seq(6), rep(1:3, each = 2)), function(y) x[y,y]))

// Lambda update, does all AFTER the E-step 
// [[Rcpp::export]]
mat lambdaUpdate(const List survtimes, const vec& ft,
				 const vec& gamma, const vec& eta, 
				 const List K, const List S,
				 const mat& b1, const mat& b2, const mat& b3,
				 const int id, 
				 const vec& w, const vec& v, const int nodes){
	mat store = zeros<mat>(ft.size(), id); // Initialise the matrix
	// Start loop over i subjects
	for(int i = 0; i < id; i++){
		vec survtimes_i = survtimes[i];    // This id's survived time indices
		// Rcpp::Rcout << "survtimes_i = " << survtimes_i << std::endl;
		// 'Unpack' the Sigma terms
		Rcpp::List Si = S[i];   
		mat S1 = as<mat>(Si[0]);
		mat S2 = as<mat>(Si[1]);
		mat S3 = as<mat>(Si[2]);  
		rowvec b1i = b1.row(i);
		// Rcpp::Rcout << "b1i = " << b1i << std::endl; 
		rowvec b2i = b2.row(i);
		rowvec b3i = b3.row(i);
		// Rcpp::Rcout << "S1 = " << S1 << std::endl;        
		rowvec Ki = K[i];                   // This id's K
		// Rcpp::Rcout << "K = " << Ki << std::endl; 
		for(int j = 0; j < survtimes_i.size(); j++){
			// Rcpp::Rcout << "ft[j] = " << ft[j] << std::endl;
			rowvec Fst = NumericVector::create(1.0, ft[j]);
			// Rcpp::Rcout << "Fst = " << Fst << std::endl;
			double tau = as_scalar(sqrt(
				pow(gamma[0], 2.0) * Fst * S1 * Fst.t() + 
				pow(gamma[1], 2.0) * Fst * S2 * Fst.t() + 
				pow(gamma[2], 2.0) * Fst * S3 * Fst.t()
			));
			// Rcpp::Rcout << "tau = " << tau << std::endl;
			// Rcpp::Rcout << "Ki * eta = " << Ki * eta << std::endl;
			// Rcpp::Rcout << "gamma[0] * b1i = " << gamma[0] * b1i << std::endl;
			// Rcpp::Rcout << "RHS = " << gamma[0] * b1i + gamma[1] * b2i + gamma[2] * b3i<< std::endl;
			double mu = as_scalar(exp(Ki * eta + Fst * (gamma[0] * b1i + 
											  gamma[1] * b2i + 
											  gamma[2] * b3i).t()));
			// Rcpp::Rcout << "mu = " << mu << std::endl;								  
			for(int k = 0; k < nodes; k++){
				store(j,i) += as_scalar(w[k] * mu * exp(v[k] * tau));
			}
		}
	}
	
	return store;
}
