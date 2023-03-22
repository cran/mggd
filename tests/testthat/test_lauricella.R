# Dimension p = 2

p <- 2
beta1 <- 1
lambda <- c(1, 2)
lauricella2_16 <- lauricella(-beta1, 0.5, p/2, 1 - lambda[1]/lambda[2], eps = 1e-16)
lauric2exact <- 1/2*(1 + lambda[1]/lambda[2])

test_that("lauricella works (dim 2)", {
  expect_equal(attr(lauricella2_16, "eps"), 1e-16)
  expect_equal(as.numeric(lauricella2_16), lauric2exact)
})

# Dimension p = 3

p <- 3
beta1 <- 1
lambda <- c(1, 1, 2)
lauricella3_16 <- lauricella(-beta1, c(0.5, 0.5), p/2, 1 - lambda[1:2]/lambda[3], eps = 1e-16)
lauric3exact <- 1/3*(1 + 2*lambda[1]/lambda[3])

test_that("lauricella works (dim 2)", {
  expect_equal(attr(lauricella2_16, "eps"), 1e-16)
  expect_equal(as.numeric(lauricella2_16), lauric2exact)
})

# Dimension p = 4

p <- 4
beta1 <- 1
lambda <- c(1, 1, 1, 2)
lauricella4_16 <- lauricella(-beta1, c(0.5, 0.5, 0.5), p/2, 1 - lambda[1:3]/lambda[4], eps = 1e-16)
lauric4exact <- 1/4*(1 + 3*lambda[1]/lambda[4])

test_that("lauricella works (dim 2)", {
  expect_equal(attr(lauricella2_16, "eps"), 1e-16)
  expect_equal(as.numeric(lauricella2_16), lauric2exact)
})
