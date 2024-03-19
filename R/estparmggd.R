estparmggd <- function(x, eps = 1e-6, display = FALSE, plot = display) {
  #' Estimation of the Parameters of a Multivariate Generalized Gaussian Distribution
  #'
  #' Estimation of the mean vector, dispersion matrix and shape parameter of a multivariate generalized Gaussian distribution (MGGD).
  #'
  #' @aliases estparmggd
  #'
  #' @usage estparmggd(x, eps = 1e-6, display = FALSE, plot = display)
  #' @param x numeric matrix or data frame.
  #' @param eps numeric. Precision for the estimation of the beta parameter.
  #' @param display logical. When \code{TRUE} the value of the \code{beta} parameter at each iteration is printed.
  #' @param plot logical. When \code{TRUE} the successive values of the \code{beta} parameter are plotted, allowing to visualise its convergence.
  #' @return A list of 3 elements:
  #' \itemize{
  #' \item \code{mu} the mean vector.
  #' \item \code{Sigma}: symmetric positive-definite matrix. The dispersion matrix.
  #' \item \code{beta} non-negative numeric value. The shape parameter.
  #' }
  #' with two attributes \code{attr(, "epsilon")} (precision of the result) and \code{attr(, "k")} (number of iterations).
  #' 
  #' @details The \eqn{\mu} parameter is the mean vector of \code{x}.
  #' 
  #' The dispersion matrix \eqn{\Sigma} and shape parameter: \eqn{\beta} are computed
  #' using the method presented in Pascal et al., using an iterative algorithm.
  #' 
  #' The precision for the estimation of \code{beta} is given by the \code{eps} parameter.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' 
  #' @references  F. Pascal, L. Bombrun, J.Y. Tourneret, Y. Berthoumieu. Parameter Estimation For Multivariate Generalized Gaussian Distribution.
  #' IEEE Trans. Signal Processing, vol. 61 no. 23, p. 5960-5971, Dec. 2013.
  #' \doi{DOI: 10.1109/TSP.2013.2282909}
  #'
  #' @seealso \code{\link{dmggd}}: probability density of a MGGD.
  #' 
  #' \code{\link{rmggd}}: random generation from a MGGD.
  #' 
  #' @examples
  #' mu <- c(0, 1, 4)
  #' Sigma <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
  #' beta <- 0.74
  #' x <- rmggd(100, mu, Sigma, beta)
  #' 
  #' # Estimation of the parameters
  #' estparmggd(x)
  #'
  #' @importFrom graphics plot
  #' @export
  
  if (is.vector(x))
    x <- cbind(x)
  if (is.data.frame(x))
    x <- as.matrix(x)
  
  # Number of variables
  p <- ncol(x)
  # Number of observations
  n <- nrow(x)
  
  # Mean vector
  mu <- colMeans(x)
  
  # Centering the data
  x <- scale(x, center = TRUE, scale = FALSE)
  
  Sigma <- diag(1, nrow = p, ncol = p)
  beta <- 0.1
  beta_1 <- Inf
  if (plot) betaseq <- numeric(0)
  invSigma <- matrix(0, nrow = n, ncol = p)
  
  if (display)
    cat(beta, "  ")
  
  k <- 0
  while (abs(beta - beta_1) > eps) {
    k <- k + 1
    sigmak <- matrix(0, nrow = p, ncol = p)
    invSigma <- x %*% solve(Sigma)
    # u <- numeric(length = n)
    u <- apply(invSigma * x, 1, sum)
    for (i in 1:n) {
      sigmak <- sigmak + u[i]^(beta-1) * t(x[i, , drop = FALSE]) %*% x[i, , drop = FALSE]
    }
    Sigma <- 1/n * sigmak
    Sigma <- (p/sum(diag(Sigma)))*Sigma
    
    beta_1 <- beta
    beta <- .estbeta(u, beta, p, eps = eps)
    if (plot) betaseq <- c(betaseq, beta)
    
    if (display)
      cat(beta, "  ")
  }
  
  if (plot)
    plot(betaseq, type = "l")
  
  if (display)
    cat("\n")
  
  m <- (beta/(p*n) * sum(u^beta))^(1/beta)
  Sigma <- m*Sigma
  
  # Returned restult
  result <- list(mu = mu, Sigma = Sigma, beta = beta)
  
  # Attributes:
  # precision of the result
  attr(result, "epsilon") <- eps
  # number of iterations
  attr(result, "k") <- k
  
  return(result)
}

.estbeta <- function(u, beta0, p, eps = .Machine$double.eps) {
  #' @importFrom stats uniroot
  
  # Search the zero of the function in ]0;beta0]
  
  N <- length(u)
  # return(
  #   pracma::fsolve(
  #     function(z) {
  #       p*N/(2*sum(u^z))*sum(log(u)*u^z) - p*N/(2*z)*(log(2) + digamma(p/(2*z))) -
  #         N - p*N/(2*z)*log(z/(p*N)*sum(u^z))
  #     }, beta0)$x  
  #      )
  result <- uniroot(
    function(z) {
      p*N/(2*sum(u^z))*sum(log(u)*u^z) - p*N/(2*z)*(log(2) + digamma(p/(2*z))) -
        N - p*N/(2*z)*log(z/(p*N)*sum(u^z))
    }, c(eps, 2*ceiling(beta0))
  )$root
  return(result)
}
