#' This function takes as an input a pxn matrix in which each column is a sampling unit
#' and each row is one marker and outputs the n-by-n quadratic kernel matrix
#' @param X is a pxn matrix in which each column is a sampling unit and each row is one marker
#' @details We obtain a quadratic kernel
#' @export quad.kernel
#' @author Arnab Maity

### Quadratic kernel
quad.kernel = function(X){  return((1+t(X)%*%X)^2)  }
