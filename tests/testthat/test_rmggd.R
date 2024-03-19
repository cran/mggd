mu1 <- c(0, 0, 0)
Sigma1 <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
beta1 <- 0.74

x1 <- rmggd(10000, mu1, Sigma1, beta1)

mu2 <- c(1, 2, 4)
Sigma2 <- matrix(c(1, 0.3, 0.2, 0.3, 0.5, 0.1, 0.2, 0.1, 0.7), nrow = 3)
beta2 <- 2.55

x2 <- rmggd(10000, mu2, Sigma2, beta2)

test_that("rmggd works",
  {
    expect_true(is.matrix(x1))
    expect_equal(dim(x1), c(10000, 3))
    expect_type(x1, "double")
    
    expect_true(is.matrix(x2))
    expect_equal(dim(x2), c(10000, 3))
    expect_type(x2, "double")
  }
)
