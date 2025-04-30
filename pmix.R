
pmix <- function(x, pi, mu, sigma) {
  pi * pnorm(x, mu[1], sigma[1]) + (1 - pi) * pnorm(x, mu[2], sigma[2])
}


