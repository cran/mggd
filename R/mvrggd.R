mvrggd <- function(n = 1, mu, Sigma, beta, tol = 1e-6) {
  #' Simulate from a Multivariate Generalized Gaussian Distribution
  #'
  #' Produces one or more samples from a multivariate (\eqn{p} variables) generalized Gaussian distribution (MGGD).
  #'
  #' @aliases mvrggd
  #'
  #' @usage mvrggd(n = 1 , mu, Sigma, beta, tol = 1e-6)
  #' @param n integer. Number of observations.
  #' @param mu length \eqn{p} numeric vector. The mean vector.
  #' @param Sigma symmetric, positive-definite square matrix of order \eqn{p}. The dispersion matrix.
  #' @param beta positive real number. The shape of the distribution.
  #' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma.
  #' @return A matrix with \eqn{p} columns and \code{n} rows.
  #' 
  #' @details A sample from a centered MGGD with dispersion matrix \eqn{\Sigma}
  #' and shape parameter \eqn{\beta} can be generated using:
  #' 
  #' \eqn{X = \tau \Sigma^{1/2} u}
  #' 
  #' where \eqn{u} is a random vector uniformly distributed on the unit sphere and
  #' \eqn{\tau} is such that \eqn{\tau^{2\beta}} is generated from a distribution
  #' \eqn{\Gamma(\frac{p}{2\beta}, 2)}.
  #' 
  #' This property is used to generate a sample from a MGGD.
  #'
  #' @author Nizar Bouhlel, Pierre Santagostini
  #' @references E. Gomez, M. Gomez-Villegas, H. Marin. A Multivariate Generalization of the Power Exponential Family of Distribution.
  #' Commun. Statist. 1998, Theory Methods, col. 27, no. 23, p 589-600.
  #' \doi{10.1080/03610929808832115}
  #' 
  #' @seealso \code{\link{mvdggd}}: probability density of a MGGD..
  #' 
  #' \code{\link{estparmvggd}}: estimation of the parameters of a MGGD.
  #'
  #' @examples
  #' mu <- c(0, 0, 0)
  #' Sigma <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
  #' beta <- 0.74
  #' mvrggd(100, mu, Sigma, beta)
  #'
  #' @importFrom MASS mvrnorm
  #' @importFrom stats rgamma
  #' @export

  # Number of variables
  p <- length(mu)

  # Sigma1 and Sigma2 must be square matrices with p rows and p columns
  if (nrow(Sigma) != p | ncol(Sigma) != p)
    stop("Sigma must be a square matrix with size equal to length(mu).")

  # IS Sigma symmetric?
  if (!isSymmetric(Sigma))
    stop("Sigma must be a symmetric, positive-definite matrix.")

  # Eigenvalues and eigenvectors of Sigma
  eig <- eigen(Sigma, symmetric = TRUE)
  lambda <- eig$values
  v <- eig$vectors

  # Square root of matrix Sigma
  A <- v %*% sqrt(diag(lambda)) %*% t(v)

  # Inverse of matrix Sigma
  invSigma <- solve(Sigma)

  # Is Sigma positive-definite?
  if (any(lambda < tol * max(abs(lambda))))
    stop("Sigma must be a symmetric, positive-definite matrix.")

  # Is beta non-negative?
  if (beta < .Machine$double.eps)
    stop("beta must be positive.")

  # r^(2*beta) is generated from a Gamma distribution with shape p/(2*beta) and scale = 2
  r <- rgamma(n, shape = p/(2*beta), scale = 2)^(1/(2*beta))

  # n samples from a multivariate Gaussian distribution with mean 0
  # and covariance matrix = identity matrix
  U <- mvrnorm(n, mu = rep(0, p), Sigma = diag(1, nrow = p), tol = tol)
  if (n == 1) U <- rbind(U)
  # Euclidian norm of each row (observation)
  normU <- sqrt(apply(U^2, 1, sum))
  # Normalize U
  U <- U / matrix(normU, nrow = n, ncol = p, byrow = FALSE)
  
  # Sample of the MGGD with mean 0, dispersion matrix Sigma and shape parameter beta
  result <- matrix(r, nrow = n, ncol = p, byrow = FALSE) * ( U %*% t(A) )
  
  # Sample of the MGGD with mean vector mu (add mu to each observation)
  result <- result + matrix(mu, nrow = n, ncol = p, byrow = TRUE)

  return(result)
}
