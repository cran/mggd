# Dimension p = 3

beta1 <- 0.74
beta2 <- 0.55
Sigma1 <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
Sigma2 <- matrix(c(1, 0.3, 0.2, 0.3, 0.5, 0.1, 0.2, 0.1, 0.7), nrow = 3)

kl12Matlab <- c("15" = 1.662970738470729,
                "20" = 1.662958569757062,
                "25" = 1.662956887529061,
                "30" = 1.662956617045686,
                "35" = 1.662956568985150,
                "40" = 1.662956559826430)

kl12_5 <- kldiv(Sigma1, beta1, Sigma2, beta2, eps = 1e-5)
kl12_6 <- kldiv(Sigma1, beta1, Sigma2, beta2, eps = 1e-6)
kl12_7 <- kldiv(Sigma1, beta1, Sigma2, beta2, eps = 1e-7)
kl12_8 <- kldiv(Sigma1, beta1, Sigma2, beta2, eps = 1e-8)

test_that("kl works (dim 3)", {
  expect_equal(attr(kl12_5, "eps"), 1e-5)
  expect_equal(attr(kl12_6, "eps"), 1e-6)
  expect_equal(attr(kl12_7, "eps"), 1e-7)
  expect_equal(attr(kl12_8, "eps"), 1e-8)
  
  expect_equal(round(as.numeric(kl12_5), 16),
               unname(kl12Matlab[as.character(attr(kl12_5, "k"))]))
  expect_equal(round(as.numeric(kl12_6), 16),
               unname(kl12Matlab[as.character(attr(kl12_6, "k"))]))
  expect_equal(round(as.numeric(kl12_7), 16),
               unname(kl12Matlab[as.character(attr(kl12_7, "k"))]))
})

# Dimension p = 4

beta1 <- 0.74
beta2 <- 0.55
Sigma1 <- matrix(c(0.8, 0.3, 0.2, 0, 0.3, 0.2, 0.1, 0, 0.2, 0.1, 0.2, 0, 0, 0, 0, 1), nrow = 4)
Sigma2 <- matrix(c(1, 0.3, 0.2, 0, 0.3, 0.5, 0.1, 0, 0.2, 0.1, 0.7, 0, 0, 0, 0, 1), nrow = 4)
kl12Matlab <- 2.017074268695473
kl21Matlab <- 9.567479623754206

kl12_9 <- kldiv(Sigma1, beta1, Sigma2, beta2, eps = 5e-9)
kl21_9 <- kldiv(Sigma2, beta2, Sigma1, beta1, eps = 5e-9)

test_that("kl works (dim 4)", {
  expect_equal(attr(kl12_9, "eps"), 5e-9)
  expect_equal(attr(kl21_9, "eps"), 5e-9)
  
  expect_equal(as.numeric(kl12_9), kl12Matlab)
  expect_equal(as.numeric(kl21_9), kl21Matlab)
})
