#' This function performs Global testig (H0: h1 =...= hp = 0) in Multivariate kernel machine regression
#' @param Y is the pxn matrix, each row is one response, each column is one subject
#' @param X.list is a list where each item of list is a qxn parametric design matrix,
#' where each row is one covariate possibly different for all items in Y
#' @param Sigma is the pxp covariance matrix of the errors
#' @param K.list is a px1 list of kernel matrices for each response
#' @param n.sim is the number of simulation to obtain the null distribution of the test statistic
#' @return test.stat is the test statistic that we obtain
#' @return p.value is the p-value of the test
#' @details We perform Global testig (H0: h1 =...= hp = 0) in Multivariate kernel machine regression
#' @export Global_test_MKMR
#' @author Arnab Maity, Mityl Biswas

Global_test_MKMR = function(Y, X.list, Sigma, K.list, n.sim = 10000){

  n = ncol(Y) ## number of subjects
  p = ncol(Sigma) ## number of response variables per subjects

  ## Identity of size n
  I.n = diag(1, nrow=n)

  # big Y-vector
  Y.til = matrix(t(Y), ncol=1)

  ## inverse of sigma
  iSigma = solve(Sigma)

  ## square root of sigma-inverse
  sig.eig <- eigen(iSigma)
  if (length(sig.eig$values) == 1)
  {  isig.sqrt <- sig.eig$vectors %*% sqrt(sig.eig$values) %*% solve(sig.eig$vectors)
  }else
    isig.sqrt <- sig.eig$vectors %*% diag(sqrt(sig.eig$values)) %*% solve(sig.eig$vectors)


  ## Covariance of the entire Y vector
  iSigma.til = kronecker(iSigma, I.n)

  ## overall Kernel matrix
  K.mat = matrix(0, ncol=n*p, nrow=n*p)
  for (k.ind in 1:p)
  {
    ind.set = c(1:n) + n*(k.ind-1)
    K.mat[ind.set, ind.set] = K.list[[k.ind]]

  }


  ## Overall parametric design matrix (X in X-transpose-beta)
  q = nrow(X.list[[1]])
  X.mat = matrix(0, nrow=q*p, ncol=n*p)
  for (x.ind in 1:p)
  {
    ind.set2 = c(1:n) + n*(x.ind-1)
    ind.set1 = c(1:q) + q*(x.ind-1)
    X.mat[ind.set1, ind.set2] = X.list[[x.ind]]

  }


  ## Hat matrix for beta
  G = solve(X.mat %*% iSigma.til %*% t(X.mat)) %*% X.mat %*% iSigma.til

  ## estimate under H0
  beta.hat = G%*% Y.til

  ## prediction of Y (under null)
  pred.Y = t(X.mat)%*%beta.hat

  ## Test statistics
  test.stat.our = t(Y.til-pred.Y) %*%    iSigma.til%*%K.mat%*%iSigma.til      %*% (Y.til-pred.Y)

  ## Simulation based test and p-value
  P0 = iSigma.til - iSigma.til%*%t(X.mat)%*%G
  P0.eig = eigen(P0)
  vtmp = P0.eig$values * (P0.eig$values > 0)
  P0.sqrt = P0.eig$vectors %*% diag(sqrt(vtmp)) %*% solve(P0.eig$vectors)


  M.mat = P0.sqrt %*% K.mat %*% P0.sqrt
  M.eig = matrix(Re(eigen(M.mat)$values), ncol=1)

  sim.tmp = rep(NA, n.sim)

  sim.r = matrix(rchisq(n.sim*nrow(M.eig), df=1), ncol=nrow(M.eig))
  sim.tmp = sim.r %*% M.eig

  pv = mean(c(sim.tmp) > c(test.stat.our))

  return(list(test.stat = c(test.stat.our), p.value = pv))

}
