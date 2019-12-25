#' Ridge Regression
#'
#' \code{ridge_regression} returns the ridge regression coefficient estimates
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
ridge_regression <- function(y, X, lambda) {
  if((nrow(X)==nrow(as.matrix(y)))==FALSE) stop("The response variable and design matrix are non-conformable")
  if(all(lambda >= 0) == FALSE) stop("The tuning parameters are negative")
  # for(i in 1:length(lambda)){betaridge <- solve((xTx + lambda[i]*diag(ncol(X))))%*%xTy
  # Coeffmat <- cbind(Coeffmat, betaridge)
  # #Coeffmat[,i]<- betaridge
  # }
  # sapply(lambda, function(lambda){solve((xTx + lambda*diag(ncol(X))))%*%xTy}))
  Coeffmat <- matrix(NA, nrow=ncol(X), ncol=0)
  xxT <- X%*%t(X)
  for(i in 1:length(lambda)){
    betaridge <- t(X)%*%solve((xxT + lambda[i]*diag(nrow(X))))%*%y
    Coeffmat <- cbind(Coeffmat, betaridge)
  }
  Coeffmat
}

#' Generalized Cross Validation
#'
#' \code{gcv} returns the leave-one-out
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
gcv <- function(y, X, lambda) {
  if((nrow(X)==nrow(as.matrix(y)))==FALSE) stop("The response variable and design matrix are non-conformable")
  if(all(lambda >= 0) == FALSE) stop("The tuning parameters are negative")
  n <- length(y)
  gencrossval<- matrix(NA,nrow=length(lambda), ncol=1)
  svdX<-svd(X)
  dsq<-svdX$d^2
  # Dsq <-dsq*diag(ncol(X))
  # UDsq <-svdX$U%*%Dsq
  # UTy <-t(svdX$u)%*%y
  xxT <- X%*%t(X)
  for(i in 1:length(lambda)){
    dof<-sum((dsq)/(dsq + lambda[i]))
    yhat<- xxT%*%solve(xxT+ lambda[i]*diag(nrow(X)))%*%y
    gencrossval[i,] <- sum(((y-yhat)/(1-(dof/n)))^2)/n
  }
  gencrossval
}

#' Leave One Out
#'
#' \code{leave_one_out} returns the leave-one-out
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
leave_one_out <- function(y, X, lambda) {
  if((nrow(X)==nrow(as.matrix(y)))==FALSE) stop("The response variable and design matrix are non-conformable")
  if(all(lambda >= 0) == FALSE) stop("The tuning parameters are negative")
  n <- length(y)
  CVlambda <- matrix(NA, nrow = length(lambda),ncol=1)
  xxT <- X%*%t(X)
  for(i in 1:length(lambda)){
    loos <- matrix(NA, nrow = n,ncol=1)
    hatmat <- xxT%*%solve((xxT + lambda[i]*diag(nrow(X))))
    yhat <-hatmat%*%y
    for(k in 1:n){
      yk<-y[k]
      yhatk<-yhat[k,]
      hatk<- hatmat[k,k]
      loos[k,]<- ((yk-yhatk)/(1-hatk))^2
    }
    CVlambda[i,] <-sum(loos)/n
  }
  CVlambda
}
