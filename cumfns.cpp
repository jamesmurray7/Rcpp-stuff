#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector cumprodC(NumericVector x) {
  int n = x.size();
  NumericVector total(n);
  total.fill(1);

  for(int i = 1; i < n; ++i) {
    total[i] = total[i-1] * x[i];
  }
  return total;
}

//[[Rcpp::export]]
NumericVector cumsumC(NumericVector x){
  int n = x.size();
  NumericVector total(n);
  total(0) = x(0);
  for(int i = 1; i < n; i++){
    total[i] = total[i-1] + x[i];
  }
  return total;
}
