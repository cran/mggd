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
