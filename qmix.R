qmix <- function(p, pi, mu, sigma) {
  f <- function(q) pmix(q, pi, mu, sigma) - p
  sapply(p, function(p_i) uniroot(f, lower = -10, upper = 10)$root)
}


