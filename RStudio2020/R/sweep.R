#' Sweep k
#'
#' \code{sweep_k} applies the sweep operator to a symmetric matrix
#' on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export
sweep_k <- function(A,k){
  if(length(diag(A)[diag(A)!=0]) !=ncol(A)) stop("There are some zero diagonal entries")
  if(!isSymmetric(A)) stop("A is not symmetric")
  hatA<-matrix(NA, nrow=nrow(A), ncol=ncol(A))
  for (i in 1:dim(A)[1]){
    for (j in 1:dim(A)[2]){
      if((j==k) | (i==k)){
        hatA[i,k] <- A[i,k]/A[k,k]
        hatA[k,j] <- A[k,j]/A[k,k]
      }else{
        hatA[i,j] <- A[i,j]-(A[i,k]*A[k,j])/A[k,k]
      }

    }
    hatA[k,k]<- -1/A[k,k]
  }
  return(hatA)
}

#' Inverse Sweep k
#'
#' \code{isweep_k} applies the inverse sweep operator to a symmetric matrix
#' on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export
isweep_k <- function(A,k){
  if(length(diag(A)[diag(A)!=0]) !=ncol(A)) stop("There are some zero diagonal entries")
  if(!isSymmetric(A)) stop("A is not symmetric")
  hatA<-matrix(NA, nrow=nrow(A), ncol=ncol(A))
  for (i in 1:dim(A)[1]){
    for (j in 1:dim(A)[2]){
      if((j==k) | (i==k)){
        hatA[i,k] <- -A[i,k]/A[k,k]
        hatA[k,j] <- -A[k,j]/A[k,k]
      }else{
        hatA[i,j] <- A[i,j]-(A[i,k]*A[k,j])/A[k,k]
      }

    }
    hatA[k,k]<- -1/A[k,k]
  }
  return(hatA)
}


#' Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep
#' @export
sweep <- function(A, k=NULL){
  #Reassigning null k to index all
  if(is.null(k)==1){
    k<-1:ncol(A)
  }
  for(i in k){
    A = sweep_k(A,i)
  }
  return(A)
}

#' Inverse Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep
#' @export
isweep <- function(A, k=NULL){
  #Reassigning null k to index all
  if(is.null(k)==1){
    k<-1:ncol(A)
  }
  for(i in k){
    A = isweep_k(A,i)
  }
  return(A)
}
