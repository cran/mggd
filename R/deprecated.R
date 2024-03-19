#' Deprecated functions
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Use [kldggd()] instead of `kldiv()`.
#' 
#' @importFrom lifecycle deprecate_soft
#' @export
#' @keywords internal
#' @name deprecated
kldiv <- function(sigma1, beta1, sigma2, beta2, eps = 1e-06) {
  deprecate_soft("1.1.0", "kldiv()", "kldggd()")
                 
  kldggd(sigma1, beta1, sigma2, beta2, eps)
}

#' @description
#' Use [dmggd()] instead of `mvdggd()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated
mvdggd <- function(x, mu, Sigma, beta, tol = 1e-6) {
  deprecate_soft("1.2.0", "mvdggd()", "dmggd()")
  
  dmggd(x, mu, Sigma, beta, tol)
}

#' @description
#' Use [rmggd()] instead of `mvrggd()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated
mvrggd <- function(n = 1, mu, Sigma, beta, tol = 1e-6) {
  deprecate_soft("1.2.0", "mvrggd()", "rmggd()")
  
  rmggd(n, mu, Sigma, beta, tol)
}

#' @description
#' Use [estparmggd()] instead of `estparmvggd()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated
estparmvggd <- function(x, eps = 1e-6, display = FALSE, plot = display) {
  deprecate_soft("1.2.0", "estparmvggd()", "estparmggd()")
  
  estparmggd(x, eps, display, plot)
}

#' @description
#' Use [plotmggd()] instead of `plotmvggd()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated
plotmvggd <- function(mu, Sigma, beta, xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
                      ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]), n = 101,
                      xvals = NULL, yvals = NULL, xlab = "x", ylab = "y",
                      zlab = "f(x,y)", col = "gray", tol = 1e-6, ...) {
  deprecate_soft("1.2.0", "plotmvggd()", "plotmggd()")
  
  plotmggd(mu, Sigma, beta, xlim, ylim, n, xvals, yvals, xlab, ylab,
           zlab, col, tol, ...)
}

#' @description
#' Use [contourmggd()] instead of `contourmvggd()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated
contourmvggd <- function(mu, Sigma, beta,
                         xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
                         ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]),
                         zlim = NULL, npt = 30, nx = npt, ny = npt,
                         main = "Multivariate generalised Gaussian density",
                         sub = NULL, nlevels = 10,
                         levels = pretty(zlim, nlevels),
                         tol = 1e-6, ...) {
  deprecate_soft("1.2.0", "contourmvggd()", "contourmggd()")
  
  contourmggd(mu, Sigma, beta, xlim, ylim, zlim, npt, nx, ny,
              main, sub, nlevels, levels, tol = 1e-6, ...)
}
