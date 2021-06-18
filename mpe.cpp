#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double mpe(List model){
  if (!model.inherits("lm")) stop("Input must be lm object");
  
  NumericVector resid = as<NumericVector>(model["residuals"]);
  NumericVector fitted = as<NumericVector>(model["fitted.values"]);
  
  int n = resid.size();
  double err = 0;
  for(int i = 0; i < n; i++){
    err +=resid[i]/(fitted[i]+resid[i]);
  }
  return err/n;
}
