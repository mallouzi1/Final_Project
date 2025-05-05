
# Mixture of Normals Functions

dmix <- function(x, pi, mu, sigma) {
  pi[1] * dnorm(x, mean = mu[1], sd = sigma[1]) + 
  pi[2] * dnorm(x, mean = mu[2], sd = sigma[2])
}

pmix <- function(x, pi, mu, sigma) {
  pi[1] * pnorm(x, mean = mu[1], sd = sigma[1]) + 
  pi[2] * pnorm(x, mean = mu[2], sd = sigma[2])
}

qmix <- function(p, pi, mu, sigma) {
  f <- function(q) pmix(q, pi, mu, sigma) - p
  sapply(p, function(p_i) uniroot(f, lower = -10, upper = 10)$root)
}

rmix <- function(n, pi, mu, sigma) {
  z <- rbinom(n, 1, pi[1])
  rnorm(n, mean = ifelse(z == 1, mu[1], mu[2]), 
           sd = ifelse(z == 1, sigma[1], sigma[2]))
}

# Visualization and Testing

# Parameters for testing
pi <- c(0.4, 0.6)
mu <- c(-2, 2)
sigma <- c(1, 1)

# Generate sample data and plot
set.seed(123)
samples <- rmix(1000, pi, mu, sigma)
hist(samples, breaks = 30, probability = TRUE, col = 'skyblue', main = "Mixture of Normals", xlab = "x")
curve(dmix(x, pi, mu, sigma), add = TRUE, col = "red", lwd = 2)

# Plot CDF and overlay qmix quantiles
p_vals <- seq(0.01, 0.99, length.out = 100)
q_vals <- qmix(p_vals, pi, mu, sigma)
cdf_vals <- pmix(q_vals, pi, mu, sigma)

plot(q_vals, cdf_vals, type = "l", col = "blue", lwd = 2, main = "qmix vs pmix", xlab = "Quantile", ylab = "CDF")
abline(0, 1, col = "gray", lty = 2)

# Test that pmix(qmix(p)) ~ p
print(all.equal(cdf_vals, p_vals, tolerance = 1e-5))
