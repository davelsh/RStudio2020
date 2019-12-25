library(testthat)
library(pracma)

###Check should run as these should throw errors
test_that("Errors",{
  expect_error(RStudio2020::ridge_regression(c(1,2), rand(10,10),c(1,2)))
  expect_error(RStudio2020::ridge_regression(rnorm(2), rand(2,2),c(-1,1)))
})

#Check should have no errors if this works
test_that("Showing Normal Equations correct", {
  y<-as.matrix(rnorm(10))
  X<-rand(10,20)
  lambda <-1
  betas <- RStudio2020::ridge_regression(y,X,lambda)
  expect_equal((t(X)%*%X + lambda*diag(ncol(X)))%*%betas, t(X)%*%y)
})
