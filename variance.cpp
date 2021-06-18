#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}

//[[Rcpp::export]]
double varC(NumericVector x){
  double mean = meanC(x);
  double sum = 0;
  int n = x.size();
  for(int i = 0; i < n; i++){
    sum += pow(x[i]-mean, 2.0);
  }
  return sum/(n-1);
}
