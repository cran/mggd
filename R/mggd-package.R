#' Tools for Multivariate Generalized Gaussian Distributions
#'
#' This package provides tools for multivariate generalized Gaussian distributions (MGGD):
#' \itemize{
#' \item Calculation of distances/divergences between multivariate generalized Gaussian distributions:
#' \itemize{
#' \item Kullback-Leibler divergence: \code{\link{kldggd}}
#' }
#' \item Tools for MGGD:
#' \itemize{
#' \item Probability density: \code{\link{dmggd}}
#' \item Estimation of the parameters: \code{\link{estparmggd}}
#' \item Simulation from a MGGD: \code{\link{rmggd}}
#' \item Plot of the density of a MGGD with 2 variables: \code{\link{plotmggd}}, \code{\link{contourmggd}}
#' }
#' }
#'
#' @name mggd-package
#' @aliases mggd-package mggd
#' @docType package
#' @author Pierre Santagostini <pierre.santagostini@agrocampus-ouest.fr>,
#' Nizar Bouhlel <nizar.bouhlel@agrocampus-ouest.fr>
#' @references N. Bouhlel, A. Dziri, Kullback-Leibler Divergence Between Multivariate Generalized Gaussian Distributions.
#' IEEE Signal Processing Letters, vol. 26 no. 7, July 2019.
#' \doi{10.1109/LSP.2019.2915000}
#' 
#' E. Gomez, M. Gomez-Villegas, H. Marin. A Multivariate Generalization of the Power Exponential Family of Distribution.
#' Commun. Statist. 1998, Theory Methods, col. 27, no. 23, p 589-600.
#' \doi{10.1080/03610929808832115}
#'
#' F. Pascal, L. Bombrun, J.Y. Tourneret, Y. Berthoumieu. Parameter Estimation For Multivariate Generalized Gaussian Distribution.
#' IEEE Trans. Signal Processing, vol. 61 no. 23, p. 5960-5971, Dec. 2013.
#' \doi{10.1109/TSP.2013.2282909}
#' #' @keywords internal 
"_PACKAGE"

NULL
