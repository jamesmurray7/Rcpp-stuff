# Rcpp stuff --------------------------------------------------------------#
# http://adv-r.had.co.nz/Rcpp.html <- from here
library(Rcpp)
library(microbenchmark)

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

x <- matrix(rnorm(1e3),nc = 100)
microbenchmark(
  rowSums(x), rowSumsC(x)
)

# Using sourceCpp
sourceCpp("~/Documents/Rcpp-stuff/mean.cpp")


# Exercises ---------------------------------------------------------------
# Function for variance
sourceCpp("~/Documents/Rcpp-stuff/variance.cpp")
x <- runif(1e5)
microbenchmark(var(x), varC(x))

sourceCpp("~/Documents/Rcpp-stuff/cumfns.cpp")
x <- 1:10
cumprod(x)
cumprodC(x)
microbenchmark(cumprod(x), cumprodC(x))
x <- 1:1000
cumsum(x)
cumsumC(x)
microbenchmark(cumsum(x), cumsumC(x))


# different classes -------------------------------------------------------
# List input
sourceCpp("~/Documents/Rcpp-stuff/mpe.cpp")

model <- lm(mpg~wt,mtcars)
mpe(model)


# STL ---------------------------------------------------------------------
sourceCpp("~/Documents/Rcpp-stuff/sumSTL.cpp")
sum3(1:100)
sum(1:100)
sumC(1:100)

x <- 1:5000
microbenchmark(sum(x), sumC(x), sum3(x), sum4(x))

# Gibbs sampler -----------------------------------------------------------
gibbs_r <- function(N, thin) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0
  
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rgamma(1, 3, y * y + 4)
      y <- rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)))
    }
    mat[i, ] <- c(x, y)
  }
  mat
}

cppFunction('NumericMatrix gibbs_cpp(int N, int thin) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rgamma(1, 3, 1 / (y * y + 4))[0];
      y = rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)))[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }

  return(mat);
}')

microbenchmark(
  gibbs_r(100, 10),
  gibbs_cpp(100, 10)
)


# Vectorisation -----------------------------------------------------------
# First way ...
vacc1a <- function(age, female, ily) {
  p <- 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily
  p <- p * if (female) 1.25 else 0.75
  p <- max(0, p)
  p <- min(1, p)
  p
}

vacc1 <- function(age, female, ily) {
  n <- length(age)
  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- vacc1a(age[i], female[i], ily[i])
  }
  out
}

# Second way...
vacc2 <- function(age, female, ily) {
  p <- 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily
  p <- p * ifelse(female, 1.25, 0.75)
  p <- pmax(0, p)
  p <- pmin(1, p)
  p
}

# Load Rcpp version...
sourceCpp("~/Documents/Rcpp-stuff/vacc.cpp")

# simulate and check
n <- 10000
age <- rnorm(n, mean = 50, sd = 10)
female <- sample(c(T, F), n, rep = TRUE)
ily <- sample(c(T, F), n, prob = c(0.8, 0.2), rep = TRUE)

stopifnot(
  all.equal(vacc1(age, female, ily), vacc2(age, female, ily)),
  all.equal(vacc1(age, female, ily), vacc3(age, female, ily))
)

# becnhmark
microbenchmark(
  vacc1 = vacc1(age, female, ily),
  vacc2 = vacc2(age, female, ily),
  vacc3 = vacc3(age, female, ily)
)