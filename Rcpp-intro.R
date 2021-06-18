# Rcpp stuff --------------------------------------------------------------
library(Rcpp)

# Simple addition ----
cppFunction('int add(int a,  int b, int c){
            int sum = a + b + c;
            return sum;
}')

rbenchmark::benchmark(
  'rcpp' = {add(100,200,300)},
  'R' = {100 + 200 + 300},
  replications = 1000
)

# Sign function ----
signR <- function(x) {
  if (x > 0) {
    1
  } else if (x == 0) {
    0
  } else {
    -1
  }
}

cppFunction('int signC(int x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}')

rbenchmark::benchmark(
  'rcpp' = {signC(-5)},
  'R' = {signR(-5)},
  replications = 10000
)

# Scalar input, scalar output ----
sumR <- function(x){
  total <- 0
  for(i in seq_along(x)){
    total <- total + x[i]
  }
  total
}

cppFunction('double sumC(NumericVector x){
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; i++){
    total += x[i];
  }
  return total;
}')

rbenchmark::benchmark(
  'rcpp' = {sumC(1:5000)},
  'R' = {sumR(1:5000)},
  "R-builtin" = {sum(1:5000)},
  replications = 100
)

microbenchmark::microbenchmark(
  sum(1:5000), sumC(1:5000), sumR(1:5000)
)

# Vector input, vector output ----
# Finding euclidian distance
pdistR <- function(x, ys){
  sqrt((x-ys) ^ 2)
}

cppFunction('NumericVector pdistC(double x, NumericVector ys){
  int n = ys.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++){
    out[i] = sqrt(pow(ys[i] - x, 2.0));
  }
  return out;
}')

x <- 10
y <- rnorm(10000, 10, 5)

microbenchmark::microbenchmark(
  pdistR(x, y), pdistC(x, y)
)

# Matrix input, vector output ----
cppFunction('NumericVector rowSumsC(NumericMatrix x){
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  
  for(int i = 0; i < nrow; i++){
    double total = 0;
    for(int j = 0; j < ncol; j++){
      total += x(i, j);
    }
    out[i] = total;
  }  
  return out;
}')

matrix(rnorm(1e5),50000,50000)


