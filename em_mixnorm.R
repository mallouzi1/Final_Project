#' Fit EM Algorithm for Mixture of Two Normals
#'
#' @param X A numeric vector of data points.
#' @param maxit Maximum number of iterations.
#' @param tolerance Convergence tolerance for log-likelihood.
#' @return A list with estimated parameters.
#' @export
em_mixnorm <- function(X, maxit = 1000, tolerance = 1e-6) {
  n <- length(X)
  PiOld <- 0.5
  clusterguess <- rbinom(n, 1, PiOld)
  
  mu1 <- mean(X[clusterguess == 1])
  sd1 <- sd(X[clusterguess == 1])
  mu2 <- mean(X[clusterguess != 1])
  sd2 <- sd(X[clusterguess != 1])
  Pi <- mean(clusterguess)
  
  loglik <- function(X, mu1, sd1, mu2, sd2, Pi, clust) {
    -log(sd1)*sum(clust) -
      sum((X[clust == 1] - mu1)^2)/(2*sd1^2) -
      log(sd2)*sum(!clust) -
      sum((X[clust == 0] - mu2)^2)/(2*sd2^2) +
      log(Pi)*sum(clust) + log(1 - Pi)*sum(!clust)
  }
  
  oldloglik <- loglik(X, mu1, sd1, mu2, sd2, Pi, clusterguess)
  
  for (i in 1:maxit) {
    l1 <- -log(sd1) - (X - mu1)^2 / (2 * sd1^2) + log(Pi)
    l2 <- -log(sd2) - (X - mu2)^2 / (2 * sd2^2) + log(1 - Pi)
    clusterguess <- l1 > l2
    
    mu1 <- mean(X[clusterguess == 1])
    sd1 <- sd(X[clusterguess == 1])
    mu2 <- mean(X[clusterguess != 1])
    sd2 <- sd(X[clusterguess != 1])
    Pi <- mean(clusterguess)
    
    newloglik <- loglik(X, mu1, sd1, mu2, sd2, Pi, clusterguess)
    if (abs(newloglik - oldloglik) < tolerance) break
    oldloglik <- newloglik
  }
  
  return(list(mu1 = mu1, sd1 = sd1, mu2 = mu2, sd2 = sd2, Pi = Pi, clusters = clusterguess))
}
