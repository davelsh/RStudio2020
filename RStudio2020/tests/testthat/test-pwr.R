library(testthat)
library(Matrix)
test_that("Errors",{
  set.seed(12345)
  n <- 1e7
  nnz <- 1e-5*n
  ix <- sample(1:n, size = nnz, replace = FALSE)
  S <- Matrix(0, nrow=n, ncol=n, sparse=TRUE)
  S[ix] <- rnorm(nnz)
  S[1,1] <- 10
  u <- matrix(rnorm(n), ncol=1)
  v <- matrix(rnorm(n), ncol=2)
  x<-rnorm(n)
  k<-1000
  tol<-1e-6
  expect_error(RStudio2020::power_method_sparse_plus_low_rank(S,u,v,x,k,tol))
})
