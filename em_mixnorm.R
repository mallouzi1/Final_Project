em_mixnorm <- function(data, max.iter = 1000, tol = 1e-6) {
  n <- length(data)
  pi <- 0.5
  mu <- quantile(data, c(0.25, 0.75))
  sigma <- rep(sd(data), 2)
  
  loglik <- function() sum(log(dmix(data, pi, mu, sigma)))
  prev_ll <- loglik()
  
  for (i in 1:max.iter) {
    # E-step
    gamma <- pi * dnorm(data, mu[1], sigma[1]) /
      (pi * dnorm(data, mu[1], sigma[1]) + (1 - pi) * dnorm(data, mu[2], sigma[2]))
    
    # M-step
    pi <- mean(gamma)
    mu[1] <- sum(gamma * data) / sum(gamma)
    mu[2] <- sum((1 - gamma) * data) / sum(1 - gamma)
    sigma[1] <- sqrt(sum(gamma * (data - mu[1])^2) / sum(gamma))
    sigma[2] <- sqrt(sum((1 - gamma) * (data - mu[2])^2) / sum(1 - gamma))
    
    curr_ll <- loglik()
    if (abs(curr_ll - prev_ll) < tol) break
    prev_ll <- curr_ll
  }
  
  list(pi = pi, mu = mu, sigma = sigma, loglik = curr_ll, iter = i)
}