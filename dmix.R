dmix <- function(x, pi, mu, sigma) {
  pi * dnorm(x, mu[1], sigma[1]) + (1 - pi) * dnorm(x, mu[2], sigma[2])
}


