kldggd <- function(Sigma1, beta1, Sigma2, beta2, eps = 1e-06) {
  #' Kullback-Leibler Divergence between Centered Multivariate generalized Gaussian Distributions
  #'
  #' Computes the Kullback- Leibler divergence between two random variables distributed
  #' according to multivariate generalized Gaussian distributions (MGGD) with zero means.
  #'
  #' @aliases kldggd
  #'
  #' @usage kldggd(Sigma1, beta1, Sigma2, beta2, eps = 1e-06)
  #' @param Sigma1 symmetric, positive-definite matrix. The dispersion matrix of the first distribution.
  #' @param beta1 positive real number. The shape parameter of the first distribution.
  #' @param Sigma2 symmetric, positive-definite matrix. The dispersion matrix of the second distribution.
  #' @param beta2 positive real number. The shape parameter of the second distribution.
  #' @param eps numeric. Precision for the computation of the Lauricella function
  #' (see \code{\link{lauricella}}). Default: 1e-06.
  #' @return A  numeric value: the Kullback-Leibler divergence between the two distributions,
  #' with two attributes \code{attr(, "epsilon")} (precision of the result of the Lauricella function;
  #' 0 if the distributions are univariate) 
  #' and \code{attr(, "k")} (number of iterations).
  #'
  #' @details Given \eqn{\mathbf{X}_1}, a random vector of \eqn{\mathbb{R}^p} (\eqn{p > 1}) distributed according to the MGGD
  #' with parameters \eqn{(\mathbf{0}, \Sigma_1, \beta_1)}
  #' and \eqn{\mathbf{X}_2}, a random vector of \eqn{\mathbb{R}^p} distributed according to the MGGD
  #' with parameters \eqn{(\mathbf{0}, \Sigma_2, \beta_2)}.
  #'
  #' The Kullback-Leibler divergence between \eqn{X_1} and \eqn{X_2} is given by:
  #' \deqn{ \displaystyle{ KL(\mathbf{X}_1||\mathbf{X}_2) = \ln{\left(\frac{\beta_1 |\Sigma_1|^{-1/2} \Gamma\left(\frac{p}{2\beta_2}\right)}{\beta_2 |\Sigma_2|^{-1/2} \Gamma\left(\frac{p}{2\beta_1}\right)}\right)} + \frac{p}{2} \left(\frac{1}{\beta_2} - \frac{1}{\beta_1}\right) \ln{2} - \frac{p}{2\beta_2} + 2^{\frac{\beta_2}{\beta_1}-1} \frac{\Gamma{\left(\frac{\beta_2}{\beta_1} + \frac{p}{\beta_1}\right)}}{\Gamma{\left(\frac{p}{2 \beta_1}\right)}} \lambda_p^{\beta_2} } }
  #' \deqn{ \displaystyle{ \times F_D^{(p-1)}\left(-\beta_1; \underbrace{\frac{1}{2},\dots,\frac{1}{2}}_{p-1}; \frac{p}{2}; 1-\frac{\lambda_{p-1}}{\lambda_p},\dots,1-\frac{\lambda_{1}}{\lambda_p}\right) } }
  #' 
  #' where \eqn{\lambda_1 < ... < \lambda_{p-1} < \lambda_p} are the eigenvalues
  #' of the matrix \eqn{\Sigma_1 \Sigma_2^{-1}}\cr
  #' and \eqn{F_D^{(p-1)}} is the Lauricella \eqn{D}-hypergeometric Function.
  #' 
  #' This computation uses the \code{\link{lauricella}} function.
  #' 
  #' When \eqn{p = 1} (univariate case):
  #' let \eqn{X_1}, a random variable distributed according to the generalized Gaussian distribution
  #' with parameters \eqn{(0, \sigma_1, \beta_1)}
  #' and \eqn{X_2}, a random variable distributed according to the generalized Gaussian distribution
  #' with parameters \eqn{(0, \sigma_2, \beta_2)}.
  #' \deqn{ KL(X_1||X_2) = \displaystyle{ \ln{\left(\frac{\frac{\beta_1}{\sqrt{\sigma_1}} \Gamma\left(\frac{1}{2\beta_2}\right)}{\frac{\beta_2}{\sqrt{\sigma_2}} \Gamma\left(\frac{1}{2\beta_1}\right)}\right)} + \frac{1}{2} \left(\frac{1}{\beta_2} - \frac{1}{\beta_1}\right) \ln{2} - \frac{1}{2\beta_2} + 2^{\frac{\beta_2}{\beta_1}-1} \frac{\Gamma{\left(\frac{\beta_2}{\beta_1} + \frac{1}{\beta_1}\right)}}{\Gamma{\left(\frac{1}{2 \beta_1}\right)}} \left(\frac{\sigma_1}{\sigma_2}\right)^{\beta_2} } }
  #' @seealso [dmggd]: probability density of a MGGD.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
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
  #' kl12 <- kldggd(Sigma1, beta1, Sigma2, beta2)
  #' kl21 <- kldggd(Sigma2, beta2, Sigma1, beta1)
  #' print(kl12)
  #' print(kl21)
  #'
  #' # Distance (symmetrized Kullback-Leibler divergence)
  #' kldist <- as.numeric(kl12) + as.numeric(kl21)
  #' print(kldist)
  #'
  #' @export
  
  # Sigma1 and Sigma2 must be matrices
  if (is.numeric(Sigma1) & !is.matrix(Sigma1))
    Sigma1 <- as.matrix(Sigma1)
  if (is.numeric(Sigma2) & !is.matrix(Sigma2))
    Sigma2 <- as.matrix(Sigma2)
  
  # Number of variables
  p <- nrow(Sigma1)
  
  # if (p == 1)
  #   stop("kldggd computes the KL divergence between multivariate probability distributions.\nSigma1 and Sigma2 must be square matrices with p > 1 rows.")
  
  # Sigma1 and Sigma2 must be square matrices with the same size
  if (ncol(Sigma1) != p | nrow(Sigma2) != p | ncol(Sigma2) != p)
    stop("Sigma1 et Sigma2 doivent must be square matrices with rank p.")
  
  # IS Sigma1 symmetric, positive-definite?
  if (!isSymmetric(Sigma1))
    stop("Sigma1 must be a symmetric, positive-definite matrix.")
  lambda1 <- eigen(Sigma1, only.values = TRUE)$values
  if (any(lambda1 < .Machine$double.eps))
    stop("Sigma1 must be a symmetric, positive-definite matrix.")
  
  # IS Sigma2 symmetric, positive-definite?
  if (!isSymmetric(Sigma2))
    stop("Sigma2 must be a symmetric, positive-definite matrix.")
  lambda2 <- eigen(Sigma2, only.values = TRUE)$values
  if (any(lambda2 < .Machine$double.eps))
    stop("Sigma2 must be a symmetric, positive-definite matrix.")
  
  if (beta1 < .Machine$double.eps)
    stop("beta1 must be non-negative.")
  if (beta2 < .Machine$double.eps)
    stop("beta2 must be non-negative.")
  
  pb1 <- p/(2*beta1)
  pb2 <- p/(2*beta2)
  
  # Eigenvalues of Sigma1 %*% inv(Sigma2)
  lambda <- sort(eigen(Sigma1 %*% solve(Sigma2), only.values = TRUE)$values, decreasing = FALSE)
  
  if (p > 1) {
    lauric <- lauricella(a = -beta2, b = rep(0.5, p-1), g = p/2, x = 1 - lambda[(p-1):1]/lambda[p],
                         eps = eps)
    result <- as.numeric(
      log( (beta1*gamma(pb2)/sqrt(det(Sigma1))) / (beta2*gamma(pb1)/sqrt(det(Sigma2))) ) +
        p/2 * (1/beta2 - 1/beta1) * log(2) - pb1 +
        2^(beta2/beta1 - 1) * gamma(beta2/beta1 + pb1) / gamma(pb1) * lambda[p]^beta2 * lauric
    )
    attr(result, "k") <- attr(lauric, "k")
    attr(result, "epsilon") <- attr(lauric, "epsilon")
  } else {
    result <- as.numeric(
      log( (beta1*gamma(pb2)/sqrt(Sigma1)) / (beta2*gamma(pb1)/sqrt(Sigma2)) ) +
        1/2 * (1/beta2 - 1/beta1) * log(2) - pb1 +
        2^(beta2/beta1 - 1) * gamma(beta2/beta1 + pb1) / gamma(pb1) * (Sigma1/Sigma2)^beta2
    )
    attr(result, "k") <- NULL
    attr(result, "epsilon") <- 0
  }
  
  return(result)
}
