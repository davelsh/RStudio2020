#' Box Constrained Least Squares
#'
#' \code{boxls} is a wrapper function around optim for doing box constrained
#' least squares estimation.
#'
#' @param X Design matrix
#' @param y response
#' @param b initial value for regression coefficient vector
#' @param lb vector of lower bounds
#' @param ub vector of upper bounds
#' @param factr L-BFGS-B parameter controling tolerance for stopping
#' @param maxit L-BFGS-B parameter controling number of maximum iterations
#' @export
boxls <- function(X, y, b, lb, ub, factr=1e7, maxit=1e2) {
  fn<-function(b){
    norm(X%*%as.matrix(b)-y,type="2")**2
  }
  optim(b,fn=fn,
        method="L-BFGS-B",
        lower=lb,upper=ub,control = list(factr=factr, maxit=maxit))
}

#' Check Suboptimality for Box Constrained Least Squares
#'
#' \code{boxls_gap} computes the suboptimality of a regression coefficient vector.
#'
#' @param X Design matrix
#' @param y response
#' @param b putative value for regression coefficient vector
#' @param lb vector of lower bounds
#' @param ub vector of upper bounds
#' @export
boxls_gap <- function(X, y, b, lb, ub) {
  gradbeta<- t(X)%*%(X%*%b)-t(X)%*%y
  theta<-ifelse(gradbeta >0, lb, ub)
  gbeta <- t(gradbeta)%*%(b-theta)
  gbeta
}
