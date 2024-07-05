plotmggd <- function(mu, Sigma, beta, xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
                      ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]), n = 101,
                      xvals = NULL, yvals = NULL, xlab = "x", ylab = "y",
                      zlab = "f(x,y)", col = "gray", tol = 1e-6, ...) {
  #' Plot of the Bivariate Generalised Gaussian Density
  #'
  #' Plots the probability density of the generalised Gaussian distribution with 2 variables
  #' with mean vector \code{mu}, dispersion matrix \code{Sigma} and shape parameter \code{beta}.
  #'
  #' @aliases plotmggd
  #'
  #' @usage plotmggd(mu, Sigma, beta, xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
  #'                  ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]), n = 101,
  #'                  xvals = NULL, yvals = NULL, xlab = "x", ylab = "y",
  #'                  zlab = "f(x,y)", col = "gray", tol = 1e-6, ...)
  #' @param mu length 2 numeric vector.
  #' @param Sigma symmetric, positive-definite square matrix of order 2. The dispersion matrix.
  #' @param beta positive real number. The shape of the distribution.
  #' @param xlim,ylim x-and y- limits.
  #' @param n A one or two element vector giving the number of steps in the x and y grid, passed to \code{\link[rgl]{plot3d.function}}.
  #' @param xvals,yvals The values at which to evaluate \code{x} and \code{y}. If used, \code{xlim} and/or \code{ylim} are ignored.
  #' @param xlab,ylab,zlab The axis labels.
  #' @param col The color to use for the plot. See \code{\link[rgl]{plot3d.function}}.
  #' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma, for the estimation of the density. see \code{\link{dmggd}}.
  #' @param ... Additional arguments to pass to \code{\link[rgl]{plot3d.function}}.
  #' @return Returns invisibly the probability density function.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' @references E. Gomez, M. Gomez-Villegas, H. Marin. A Multivariate Generalization of the Power Exponential Family of Distribution.
  #' Commun. Statist. 1998, Theory Methods, col. 27, no. 23, p 589-600.
  #' \doi{10.1080/03610929808832115}
  #'
  #' @seealso \code{\link{contourmggd}}: contour plot of a bivariate generalised Gaussian density.
  #' 
  #' \code{\link{dmggd}}: Probability density of a multivariate generalised Gaussian distribution.
  #'
  #' @examples
  #' mu <- c(1, 4)
  #' Sigma <- matrix(c(0.8, 0.2, 0.2, 0.2), nrow = 2)
  #' beta <- 0.74
  #' plotmggd(mu, Sigma, beta)
  #'
  #' @import rgl
  #' @importFrom rgl plot3d
  #' @export
  
  if (length(mu)!=2 | nrow(Sigma)!=2 | ncol(Sigma)!=2)
    stop(paste("plotmggd only allows plotting a generalised Gaussian density with 2 variables.",
               "mu must be a length 2 numeric vector and Sigma must be a 2*2 square matrix.", sep = "\n"))
  
  # Estimation of the density
  f <- function(x) dmggd(x, mu = mu, Sigma = Sigma, beta = beta, tol = tol)
  ff <- function(x, y) sapply(1:length(x), function(i) as.numeric(f(c(x[i], y[i]))))
  
  if (length(n) == 1)
    n <- rep(n, 2)
  if (is.null(xvals))
    xvals = seq.int(min(xlim), max(xlim), length.out = n[1])
  if (is.null(yvals))
    yvals = seq.int(min(ylim), max(ylim), length.out = n[2])
  
  # Plot
  plot3d(ff, xlim = xlim, ylim = ylim, n = n, xvals = xvals, yvals = yvals,
         xlab = xlab, ylab = ylab, zlab = zlab, col = col, ...)
  
  return(invisible(f))
}
