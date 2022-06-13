kldiv <- function(sigma1, beta1, sigma2, beta2, eps = 1e-06) {
  #' Kullback-Leibler Divergence between Centered Multivariate generalized Gaussian Distributions
  #'
  #' Computes the Kullback- Leibler divergence between two random variables distributed
  #' according to multivariate generalized Gaussian distributions (MGGD) with zero means.
  #'
  #' @aliases kldiv
  #'
  #' @usage kldiv(sigma1, beta1, sigma2, beta2, eps = 1e-06)
  #' @param sigma1 symmetric, positive-definite matrix. The dispersion matrix of the first distribution.
  #' @param beta1 positive real number. The shape parameter of the first distribution.
  #' @param sigma2 symmetric, positive-definite matrix. The dispersion matrix of the second distribution.
  #' @param beta2 positive real number. The shape parameter of the second distribution.
  #' @param eps numeric. Precision for the computation of the Lauricella function
  #' (see \code{\link{lauricella}}). Default: 1e-06.
  #' @return A  numeric value: the Kullback-Leibler divergence between the two distributions,
  #' with two attributes \code{attr(, "epsilon")} (precision of the result of the Lauricella function) 
  #' and \code{attr(, "k")} (number of iterations).
  #'
  #' @details Given \eqn{X_1}, a random vector of \eqn{R^p} distributed according to the MGGD
  #' with parameters \eqn{(0, \Sigma_1, \beta_1)}
  #' and \eqn{X_2}, a random vector of \eqn{R^p} distributed according to the MGGD
  #' with parameters \eqn{(0, \Sigma_2, \beta_2)}.
  #'
  #' The Kullback-Leibler divergence between \eqn{X_1} and \eqn{X_2} is given by:
  #' 
  #' \eqn{ KL(X_1||X_2) = ln{\left(\frac{\beta_1 |\Sigma_1|^{-1/2} \Gamma\left(\frac{p}{2\beta_2}\right)}{\beta_2 |\Sigma_2|^{-1/2} \Gamma\left(\frac{p}{2\beta_1}\right)}\right)} + \frac{p}{2} \left(\frac{1}{\beta_2} - \frac{1}{\beta_1}\right) ln{2} - \frac{p}{2\beta_2} + 2^{\frac{\beta_2}{\beta_1}-1} \frac{\Gamma{\left(\frac{\beta_2}{\beta_1} + \frac{p}{\beta_1}\right)}}{\Gamma{\left(\frac{p}{2 \beta_1}\right)}} \lambda_p^{\beta_2} }
  #' 
  #' \eqn{\times F_D^{(p-1)}\left(-\beta_1; \frac{1}{2},\dots,\frac{1}{2}; \frac{p}{2}; 1-\frac{\lambda_{p-1}}{\lambda_p},\dots,1-\frac{\lambda_{1}}{\lambda_p}\right) }
  #' 
  #' where \eqn{\lambda_1 < ... < \lambda_{p-1} < \lambda_p} are the eigenvalues
  #' of the matrix \eqn{\Sigma_1 \Sigma_2^{-1}}.
  #' 
  #' This computation uses the \code{\link{lauricella}} function.
  #'
  #' @author Nizar Bouhlel, Pierre Santagostini
  #' @references N. Bouhlel, A. Dziri, Kullback-Leibler Divergence Between Multivariate Generalized Gaussian Distributions.
  #' IEEE Signal Processing Letters, vol. 26 no. 7, July 2019.
  #' \doi{10.1109/LSP.2019.2915000}
  #'
  #' @examples
  #' beta1 <- 0.74
  #' beta2 <- 0.55
  #' Sigma1 <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
  #' Sigma2 <- matrix(c(1, 0.3, 0.2, 0.3, 0.5, 0.1, 0.2, 0.1, 0.7), nrow = 3)
  #'
  #' # Kullback-Leibler divergence
  #' kl12 <- kldiv(Sigma1, beta1, Sigma2, beta2)
  #' kl21 <- kldiv(Sigma2, beta2, Sigma1, beta1)
  #' print(kl12)
  #' print(kl21)
  #'
  #' # Distance (symmetrized Kullback-Leibler divergence)
  #' kldist <- as.numeric(kl12) + as.numeric(kl21)
  #' print(kldist)
  #'
  #' @export

  # Number of variables
  p <- nrow(sigma1)

  # sigma1 and sigma2 must be square matrices with the same size
  if (ncol(sigma1) != p | nrow(sigma2) != p | ncol(sigma2) != p)
    stop("sigma1 et sigma2 doivent must be square matrices with rank p.")

  # IS sigma1 symmetric, positive-definite?
  if (!isSymmetric(sigma1))
    stop("sigma1 must be a symmetric, positive-definite matrix.")
  lambda1 <- eigen(sigma1, only.values = TRUE)$values
  if (any(lambda1 < .Machine$double.eps))
    stop("sigma1 must be a symmetric, positive-definite matrix.")

  # IS sigma2 symmetric, positive-definite?
  if (!isSymmetric(sigma2))
    stop("sigma2 must be a symmetric, positive-definite matrix.")
  lambda2 <- eigen(sigma2, only.values = TRUE)$values
  if (any(lambda2 < .Machine$double.eps))
    stop("sigma2 must be a symmetric, positive-definite matrix.")

  if (beta1 < .Machine$double.eps)
    stop("beta1 must be non-negative.")
  if (beta2 < .Machine$double.eps)
    stop("beta2 must be non-negative.")

  pb1 <- p/(2*beta1)
  pb2 <- p/(2*beta2)

  # Eigenvalues of sigma1 %*% inv(sigma2)
  lambda <- sort(eigen(sigma1 %*% solve(sigma2), only.values = TRUE)$values, decreasing = FALSE)

  lauric <- lauricella(a = -beta2, b = rep(0.5, p-1), g = p/2, x = 1 - lambda[(p-1):1]/lambda[p],
                        eps = eps)
  
  # cat("lauric = ", lauric, "\n")

  result <- log( (beta1*gamma(pb2)/sqrt(det(sigma1))) / (beta2*gamma(pb1)/sqrt(det(sigma2))) ) +
      p/2 * (1/beta2 - 1/beta1) * log(2) - pb1 +
      2^(beta2/beta1 - 1) * gamma(beta2/beta1 + pb1) / gamma(pb1) * lambda[p]^beta2 * lauric
  attr(result, "k") <- attr(lauric, "k")
  attr(result, "epsilon") <- attr(lauric, "epsilon")

  return(result)
}
