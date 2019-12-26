#' Estimate auc
#'
#' \code{estimate_auc} estimates the area under the curve via Monte Carlo
#'
#' @param n number of darts to throw
#' @param a left limit
#' @param b right limit
#' @export
estimate_auc <- function(n, a, b) {
  ((b-a)/(n*sqrt(2*pi)))*sum(exp(-(runif(n, a, b))^2/2)/sqrt(2*pi)
                             >= runif(n, 0, 1/sqrt(2*pi)))
}



#' Estimate auc
#'
#' \code{estimate_auc_tol} estimates the area under the curve via Monte Carlo
#'
#' @param a left limit
#' @param b right limit
#' @param tol tolerance
#' @param p probability
#' @export
estimate_auc_tol <- function(a, b, tol, p) {
  n <- -log(p/2)/(2*(tol)^2)
  auctol <- ((b-a)/(n*sqrt(2*pi)))*sum(exp(-(runif(n, a, b))^2/2)/sqrt(2*pi)
                                       >= runif(n, 0, 1/sqrt(2*pi)))
  #Returning the estimated AUC with a specified tolerance, as well as the approximate n
  return(c(auctol, ceiling(n)))
}

#' Estimate auc modified
#'
#' \code{estimate_auc_modified} estimates the area under the curve via Monte Carlo
#'
#' @param n number of darts to throw
#' @param a left limit
#' @param b right limit
#' @param c canvas dimension variable
#' @export
estimate_auc_modified <- function(n, a, b, c=1/sqrt(2*pi)) {
  if(c< 1/sqrt(2*pi)) stop("Canvas height too low")

  ((b-a)/(n*sqrt(2*pi)))*sum(exp(-(runif(n, a, b))^2/2)/sqrt(2*pi)
                             >= runif(n, 0, c))

}
