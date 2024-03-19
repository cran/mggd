set.seed(42)

mu1 <- c(0, 0, 0)
Sigma1 <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
beta1 <- 0.74

x1 <- rmggd(10000, mu1, Sigma1, beta1)

set.seed(42)

mu2 <- c(1, 2, 4)
Sigma2 <- matrix(c(1, 0.3, 0.2, 0.3, 0.5, 0.1, 0.2, 0.1, 0.7), nrow = 3)
beta2 <- 2.55

x2 <- rmggd(10000, mu2, Sigma2, beta2)

# Esimation des paramètres
x1estim <- estparmggd(x1, display = FALSE, plot = FALSE)
x2estim <- estparmggd(x2, display = FALSE, plot = FALSE)

test_that("estparmggd works",
  {
    expect_equal(mu1, x1estim$mu, tolerance = 0.1)
    expect_equal(Sigma1, x1estim$Sigma, tolerance = 0.1)
    expect_equal(beta1, x1estim$beta, tolerance = 0.1)
    
    expect_equal(mu2, x2estim$mu, tolerance = 0.1)
    expect_equal(Sigma2, x2estim$Sigma, tolerance = 0.1)
    expect_equal(beta2, x2estim$beta, tolerance = 0.1)
  }
)
