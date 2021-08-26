#' This function performs classical parametric multivariate regression, Y = X\beta +  e
#' and estimates \beta, a covariance matrix for estimated \beta, and a vector of residuals
#' @param Y is the pxn matrix, each row is one response, each column is one subject
#' @param X.list is a list where each item of list is a qxn parametric design matrix,
#' where each row is one covariate possibly different for all items in Y
#' @param Sigma is the pxp covariance matrix of the errors
#' @return beta.hat is the estimate of beta, a qp by 1 vector, first q elements correspond to first response variable and so on
#' @return beta.cov is the estimated covariance matrix for estimated beta parameters
#' @return resid is the vector of residuals
#' @details We perform classical parametric multivariate regression, Y = X\beta +  e
#' @export Parametric_MR
#' @author Arnab Maity

###########################
###
###########################
Parametric_MR = function(Y, X.list, Sigma){

  n = ncol(Y)  ## number of subjects
  p = ncol(Sigma) ## number of response variables per subject


  ## Identity of size n
  I.n = diag(1, nrow=n)

  # big Y-vector
  Y.til = matrix(t(Y), ncol=1)

  ## inverse of sigma
  iSigma = solve(Sigma)

  ## square root of sigma-inverse
  sig.eig <- eigen(iSigma)
  isig.sqrt <- sig.eig$vectors %*% diag(sqrt(sig.eig$values)) %*% solve(sig.eig$vectors)


  ## Covariance of the entire Y vector
  Sigma.til = kronecker(Sigma, I.n)

  ## temporary var (needed later)
  iPs = kronecker(iSigma, I.n)

  ## Overall parametric design matrix (X in X-transpose-beta)
  q = nrow(X.list[[1]])
  X.mat = matrix(0, nrow=q*p, ncol=n*p)
  for (x.ind in 1:p)
  {
    ind.set2 = c(1:n) + n*(x.ind-1)
    ind.set1 = c(1:q) + q*(x.ind-1)
    X.mat[ind.set1, ind.set2] = X.list[[x.ind]]

  }

  ## hat matrix for beta
  G = solve(X.mat %*% iPs %*% t(X.mat)) %*% X.mat %*% iPs

  ## Estimates
  beta.hat = G%*% Y.til
  beta.sd = G %*% Sigma.til %*% t(G)

  ## predictions
  pred.Y = t(X.mat)%*%beta.hat

  #return
  return(list(beta.hat = beta.hat, beta.cov = beta.sd, resid = Y.til-pred.Y))

}
