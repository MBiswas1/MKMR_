#' This function gives us p-values for testing if a scalar predictor is correlated to a functional response.
#' @param data is a data frame comprising: id, the individual's identifier
#'                                          y, the response of the individual
#'                                          x, the argument of the response function
#' @param Z is a matrix whose each column corresponds to a covariate of interest and
#' each row corresponds to each individual
#' @param Cov is a matrix whose each column corresponds to a nuisance covariate and
#' each row corresponds to each individual
#' @details A functional response variable Y is regressed on scalar variables X and Z to determine if X is a significantly associated with Y.
#' We obtain the three outputs from SITAR and the scores obtained from fpca with 95% pve
#' and use multivariate kernel machine regression on them to get p-values using
#' linear and quadratic kernels
#' @return pv.fin is a list of p-values obtained by using the linear
#' and quadratic kernels and applying Bonferroni's correction and Sidak's correction
#' to them by using both SITAR and fpca. A value of -1 indicates the corresponding
#' method failed to converge.
#' @export pv
#' @author Mityl Biswas

pv <- function(data, Z, Cov)
{
  require(sitar)
  require(refund)
  pv <- NULL
  pv.fin <- NULL
  res <- NULL
  res2 <- NULL
  X <- t(as.matrix(Z)) # Covariate
  K1 <- linear.kernel(X)
  K2 <- quad.kernel(X)
  K.list0 <- list(K1, K2)

  #sitar
  try({test <- sitar(x, y, id, data, df = 5, method = "REML",
                     control = nlmeControl(msMaxIter = 1000, returnObject = TRUE))
  Y <- test$coefficients$random$id[,1:3]
  out <- lm(Y ~ Cov - 1) # Nuisance
  Sigma <- cov(out$residuals) # Get Sigma
  # Test
  p <- ncol(Y)
  for(j in 1:2)
    try({
      Klist <- rep(list(K.list0[[j]]), p)
      Xlist <- rep(list(t(Cov)), p)
      res[[j]] <- Global_test_MKMR(Y = t(Y),
                                   X.list = Xlist, # list of design mats
                                   Sigma = Sigma,
                                   K.list = Klist,
                                   n.sim = 10000)$p.value})
  }  )

  #fpca
  data.pca <- data.frame(.id = data$id,
                         .index = data$x,
                         .value = data$y)
  try(
    {test.pca <- fpca.sc(ydata = data.pca, pve = 0.95)
    Y2 <- test.pca$scores
    n <- ncol(Y2)
    out2 <- lm(Y2 ~ Cov - 1) # Nuisance
    Sigma2 <- cov(as.matrix(out2$residuals)) # Get Sigma
    p2 <- ncol(Y2)
    # Test
    for(j in 1:2)
    {
      Klist2 <- rep(list(K.list0[[j]]), p2)
      Xlist2 <- rep(list(t(Cov)), p2)
      res2[[j]] <- Global_test_MKMR(Y = t(Y2),
                                    X.list = Xlist2,
                                    Sigma = Sigma2,
                                    K.list = Klist2,
                                    n.sim = 10000)$p.value
    }
    })

  pv$l.sitar <- res[[1]]
  pv$q.sitar <- res[[2]]
  pv$sit.bon <- min(unlist(res) * 2, 1)
  pv$sit.sid <- min(1-(1-unlist(res))^2)

  pv$l.fpca <- res2[[1]]
  pv$q.fpca <- res2[[2]]
  pv$pca.bon <- min(unlist(res2) * 2, 1)
  pv$pca.sid <- min(1-(1-unlist(res2))^2)

  if(is.na(pv$l.sitar) || is.null(pv$l.sitar))
  {pv.fin$l.sitar <- -1} else
    pv.fin$l.sitar <- pv$l.sitar

  if(is.na(pv$q.sitar) || is.null(pv$q.sitar))
  {pv.fin$q.sitar <- -1} else
    pv.fin$q.sitar <- pv$q.sitar

  if(is.na(pv$l.fpca) || is.null(pv$l.fpca))
  {pv.fin$l.fpca <- -1} else
    pv.fin$l.fpca <- pv$l.fpca

  if(is.na(pv$q.fpca) || is.null(pv$q.fpca))
  {pv.fin$q.fpca <- -1} else
    pv.fin$q.fpca <- pv$q.fpca

  pv$sit.bon <- min(unlist(res) * 2, 1)
  pv$sit.sid <- min(1-(1-unlist(res))^2)

  pv.fin$sit.bon <- pv$sit.bon
  pv.fin$sit.sid <- pv$sit.sid

  pv.fin$pca.bon <- pv$pca.bon
  pv.fin$pca.sid <- pv$pca.sid
  return(pv.fin)
}
