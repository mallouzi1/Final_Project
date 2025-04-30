rmix <- function(n, pi, mu, sigma) {
  z <- rbinom(n, 1, pi)
  rnorm(n, mean = ifelse(z == 1, mu[1], mu[2]), sd = ifelse(z == 1, sigma[1], sigma[2]))
}


