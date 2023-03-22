mvdggd <- function(x, mu, Sigma, beta, tol = 1e-6) {
  #' Density of a Multivariate Generalized Gaussian Distribution
  #'
  #' Density of the multivariate (\eqn{p} variables) generalized Gaussian distribution (MGGD)
  #' with mean vector \code{mu}, dispersion matrix \code{Sigma} and shape parameter \code{beta}.
  #'
  #' @aliases mvdggd
  #'
  #' @usage mvdggd(x, mu, Sigma, beta, tol = 1e-6)
  #' @param x length \eqn{p} numeric vector.
  #' @param mu length \eqn{p} numeric vector. The mean vector.
  #' @param Sigma symmetric, positive-definite square matrix of order \eqn{p}. The dispersion matrix.
  #' @param beta positive real number. The shape of the distribution.
  #' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma.
  #' @return The value of the density.
  #' 
  #' @details The density function of a multivariate generalized Gaussian distribution is given by:
  #' \deqn{ \displaystyle{ f(\mathbf{x}|\boldsymbol{\mu}, \Sigma, \beta) = \frac{\Gamma\left(\frac{p}{2}\right)}{\pi^\frac{p}{2} \Gamma\left(\frac{p}{2 \beta}\right) 2^\frac{p}{2\beta}} \frac{\beta}{|\Sigma|^\frac{1}{2}} e^{-\frac{1}{2}\left((\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x}-\boldsymbol{\mu})\right)^\beta} } }
  #' 
  #' When \eqn{p=1} (univariate case) it becomes:
  #' \deqn{ \displaystyle{ f(x|\mu, \sigma, \beta) = \frac{\Gamma\left(\frac{1}{2}\right)}{\pi^\frac{1}{2} \Gamma\left(\frac{1}{2 \beta}\right) 2^\frac{1}{2\beta}} \frac{\beta}{\sigma^\frac{1}{2}} \ e^{-\left(\frac{(x - \mu)^2}{2 \sigma}\right)^\beta} = \frac{\beta}{\Gamma\left(\frac{1}{2 \beta}\right) 2^\frac{1}{2 \beta} \sqrt{\sigma}} \ e^{-\left(\frac{(x - \mu)^2}{2 \sigma}\right)^\beta} } }
  #' 
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' @references E. Gomez, M. Gomez-Villegas, H. Marin. A Multivariate Generalization of the Power Exponential Family of Distribution.
  #' Commun. Statist. 1998, Theory Methods, col. 27, no. 23, p 589-600.
  #' \doi{10.1080/03610929808832115}
  #'
  #' @seealso \code{\link{mvrggd}}: random generation from a MGGD.
  #' 
  #' \code{\link{estparmvggd}}: estimation of the parameters of a MGGD.
  #' 
  #' \code{\link{plotmvggd}}, \code{\link{contourmvggd}}: plot of the probability density of a bivariate generalised Gaussian distribution.
  #'
  #' @examples
  #' mu <- c(0, 1, 4)
  #' Sigma <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
  #' beta <- 0.74
  #' mvdggd(c(0, 1, 4), mu, Sigma, beta)
  #' mvdggd(c(1, 2, 3), mu, Sigma, beta)
  #'
  #' @export

  # Number of variables
  p <- length(mu)
  
  # Sigma must be a matrix
  if (is.numeric(Sigma) & !is.matrix(Sigma))
    Sigma <- as.matrix(Sigma)
  
  # x must have the same length as mu
  if (length(x) != p)
    stop(paste("x does not have", p, "elements.\n x and mu must have the same length."))

  # Sigma1 and Sigma2 must be square matrices with p rows and p columns
  if (nrow(Sigma) != p | ncol(Sigma) != p)
    stop("Sigma must be a square matrix with size equal to length(mu).")

  # IS Sigma symmetric?
  if (!isSymmetric(Sigma))
    stop("Sigma must be a symmetric, positive-definite matrix.")

  # Eigenvalues and eigenvectors of Sigma
  eig <- eigen(Sigma, symmetric = TRUE)
  lambda <- eig$values

  # Is Sigma positive-definite?
  if (any(lambda < tol * max(abs(lambda))))
    stop("Sigma must be a symmetric, positive-definite matrix.")
  
  # Is beta positive?
  if (beta < .Machine$double.eps)
    stop("beta must be positive.")

  # Inverse of matrix Sigma
  invSigma <- solve(Sigma)
  
  xcent <- cbind(x - mu)
  p2b <- p/(2*beta)

  # Computation of the density
  result <- gamma(p/2) / ( pi^(p/2)*gamma(p2b)*2^p2b )
  result <- result * beta/sqrt(det(Sigma))
  result <- result * exp( -0.5*abs( t(xcent) %*% invSigma %*% xcent )^beta )
  
  return(as.numeric(result))
}
