#' Power Method for dense matrices
#'
#' \code{power_method_dense} applies the power method to estimate
#' the right singular vector of a dense matrix associated with its largest singular value.
#'
#' @param A The input matrix
#' @param x initial guess
#' @param max_iter maximum number of iterations
#' @param tol relative error stopping criterion
#' @export
power_method_dense <- function(A, x, max_iter, tol) {


  x_last<-x
  for (i in 2:max_iter){

    num<-t(A)%*%(A%*%x_last)
    x_new <- num/norm(num,type="2")
    if (norm(x_new - x_last, type="2") <= tol*norm(x_last, type="2")){break}
    x_last = x_new
  }
  x_new
}


#' Check sparsity
#'
#' \code{is.sparseMatrix} checks to ensure a matrix is a sparse matrix
#' created using sparseMatrix.
#'
#' @param x The input matrix
#' @export
is.sparseMatrix <- function(x){
  is(x, 'sparseMatrix')
}

#' Power Method for sparse matrices
#'
#' \code{power_method_sparse} applies the power method to estimate
#' the right singular vector of a sparse matrix associated with its largest singular value.
#'
#' @param A The input matrix
#' @param x initial guess
#' @param max_iter maximum number of iterations
#' @param tol relative error stopping criterion
#' @export
power_method_sparse <- function(A, x, max_iter, tol) {
  if(dbelsheiST758::is.sparseMatrix(A) == FALSE) stop("Matrix isn't sparse")
  A<-as.matrix(A)
  x_last<-x
  for (i in 2:max_iter){

    num<-t(A)%*%(A%*%x_last)
    x_new <- num/norm(num,type="2")
    if (norm(x_new - x_last, type="2") <= tol*norm(x_last, type="2")){break}
    x_last = x_new
  }
  x_new
}

#' Power Method for low rank matrices
#'
#' \code{power_method_low_rank} applies the power method to estimate
#' the right singular vector of a low rank matrix associated with its largest singular value.
#'
#' @param U The left input factor matrix
#' @param V The right input factor matrix
#' @param x initial guess
#' @param max_iter maximum number of iterations
#' @param tol relative error stopping criterion
#' @export
power_method_low_rank <- function(U, V, x, max_iter, tol) {
  if(ncol(U)!=ncol(V)) stop("U, V don't have the same number of columns" )
  A <- U%*%t(V)

  x_last<-x
  for (i in 2:max_iter){

    num<-t(A)%*%(A%*%x_last)
    x_new <- num/norm(num,type="2")
    if (norm(x_new - x_last, type="2") <= tol*norm(x_last, type="2")){break}
    x_last = x_new
  }
  x_new
}

#' Power Method for sparse + low rank matrices
#'
#' \code{power_method_sparse_plus_low_rank} applies the power method to estimate
#' the right singular vector of a sparse + low rank matrix associated with its largest singular value.
#'
#' @param S sparse input matrix term
#' @param U The left input factor matrix term
#' @param V The right input factor matrix term
#' @param x initial guess
#' @param max_iter maximum number of iterations
#' @param tol relative error stopping criterion
#' @export
power_method_sparse_plus_low_rank <- function(S, U, V, x, max_iter, tol) {
  if(ncol(U)!=ncol(V)) stop("U, V don't have the same number of columns" )
  if(dbelsheiST758::is.sparseMatrix(S) == FALSE) stop("Matrix isn't sparse")


  x_last<-x
  for (i in 2:max_iter){

    num<-as.matrix((Matrix::t(S)%*%S)%*%x_last +
                     U%*%(t(U)%*%U)%*%(t(U)%*%x_last) +
                     U%*%(t(U)%*%(S%*%x_last)) +
                     (Matrix::t(S)%*%U)%*%(t(U)%*%x_last))
    x_new <- num/norm(num,type="2")
    if (norm(x_new - x_last, type="2") <= tol*norm(x_last, type="2")){break}
    x_last = x_new
  }
  x_new
}
