#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
bool s_v_d(mat& A){
  mat U;
  vec s;
  mat V;
  return svd(U,s,V,A);
}

// [[Rcpp::export]]
List baseSVD(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "standard");
  return List::create(Named("U") = U, _["V"] = V, _["d"] = S);
}

// [[Rcpp::export]]
List dcSVD(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "dc");
  return List::create(Named("U") = U, _["V"] = V, _["d"] = S);
}

// [[Rcpp::export]]
vec sqrt_diag(const mat & X){
  List a = dcSVD(X);
  mat U = a["U"], V = a["V"];
  vec d = a["d"];
  mat D = diagmat(d);
  
  return diagvec(U * D.i() * V.t());
}

