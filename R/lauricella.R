lauricella <- function(a, b, g, x, eps = 1e-06) {
  #' Lauricella \eqn{D}-Hypergeometric Function
  #'
  #' Computes the Lauricella \eqn{D}-hypergeometric Function function.
  #'
  #' @aliases lauricella
  #'
  #' @usage lauricella(a, b, g, x, eps = 1e-06)
  #' @param a numeric.
  #' @param b numeric vector.
  #' @param g numeric.
  #' @param x numeric vector. \code{x} must have the same length as \code{b}.
  #' @param eps numeric. Precision for the nested sums (default 1e-06).
  #' @return A numeric value: the value of the Lauricella function,
  #' with two attributes \code{attr(, "epsilon")} (precision of the result) and \code{attr(, "k")} (number of iterations).
  #'
  #' @details If \eqn{n} is the length of the \eqn{b} and \code{x} vectors,
  #' the Lauricella \eqn{D}-hypergeometric Function function is given by:
  #' \deqn{\displaystyle{F_D^{(n)}\left(a, b_1, ..., b_n, g; x_1, ..., x_n\right) = \sum_{m_1 \geq 0} ... \sum_{m_n \geq 0}{ \frac{ (a)_{m_1+...+m_n}(b_1)_{m_1} ... (b_n)_{m_n} }{ (g)_{m_1+...+m_n} } \frac{x_1^{m_1}}{m_1!} ... \frac{x_n^{m_n}}{m_n!} } }}
  #' 
  #' where \eqn{(x)_p} is the Pochhammer symbol (see \code{\link{pochhammer}}).
  #' 
  #' If \eqn{|x_i| < 1, i = 1, \dots, n}, this sum converges.
  #' Otherwise there is an error.
  #' 
  #' The \code{eps} argument gives the required precision for its computation.
  #' It is the \code{attr(, "epsilon")} attribute of the returned value.
  #' 
  #' Sometimes, the convergence is too slow and the required precision cannot be reached.
  #' If this happens, the \code{attr(, "epsilon")} attribute is the precision that was really reached.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' @references N. Bouhlel, A. Dziri, Kullback-Leibler Divergence Between Multivariate Generalized Gaussian Distributions.
  #' IEEE Signal Processing Letters, vol. 26 no. 7, July 2019.
  #' \doi{10.1109/LSP.2019.2915000}
  #' @importFrom utils combn
  #' @export

  # Number of variables
  n <- length(x)
  
  # Do x and b have the same length?
  if (length(b) != n)
    stop("x and b must have the same length")
  
  # Condition for the convergence: are all abs(x) < 1 ?
  if (any(abs(x) >= 1))
    stop("The series does not converge for these x values.")
  
  access <- function(ind, tab) {
    sapply(1:n, function(ii) tab[ind[ii] + 1, ii])
  }
  
  k <- 5
  
  # M: data.frame of the indices for the nested sums
  # (i.e. all arrangements of n elements from {0:k})
  M <- expand.grid(rep(list(0:k), n))

  # Sum of the indices
  Msum <- rowSums(M)
  
  Munique <- 0:max(M)
  Msumunique <- 0:max(Msum)
  
  # Product x^{m_1} * ... * x^{m_n} for m_1 = 0...k, ...,  m_n = 0...k
  # xfact <- apply(M, 1, function(Mi) prod( x^Mi ))
  # lnfact <- apply(M, 1, function(Mi) sum(sapply(Mi, lnfactorial)))
  xfact <- as.data.frame(
    matrix(nrow = max(M) + 1, ncol = n, dimnames = list(Munique, 1:n)))
  for (i in Munique) for (j in 1:n) {
    # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
    xfact[i+1, j] <- x[j]^i
  }
  # prodxfact <- function(ind) {
  #   # prod(mapply(function(i, j) xfact[i, j], ind+1, 1:n))
  #   prod(access(ind, xfact))
  # }
  gridxfact <- expand.grid(xfact)
  prodxfact <- apply(gridxfact, 1, prod)
  
  # Logarithm of the product m_1! * ... * m_n! for m_1 = 0...k, ...,  m_n = 0...k
  # i.e. \sum_{i=0}^n{\log{m_i!}}
  # lnfact <- as.data.frame(
  #   matrix(sapply(Munique, lnfactorial), nrow = length(Munique),
  #          ncol = n, byrow = FALSE, dimnames = list(Munique, 1:n)))
  lnfact <- as.data.frame(
    matrix(lfactorial(Munique), nrow = length(Munique),
           ncol = n, byrow = FALSE, dimnames = list(Munique, 1:n)))
  # sumlnfact <- function(ind) {
  #   # sum(mapply(function(i, j) lnfact[i, j], ind+1, 1:n))
  #   sum(access(ind, lnfact))
  # }
  # Table of all combinations of the log(m_i)
  gridlnfact <- expand.grid(lnfact)
  # Sum of the logarithms
  sumlnfact <- rowSums(gridlnfact)
  
  # Logarithms of pochhammer(a, m_1+...+m_n) for m_1 = 0...k, ...,  m_n = 0...k
  # lnapoch <- sapply(Msum, function(j) lnpochhammer(a, j))
  lnapoch <- sapply(Msumunique, function(j) lnpochhammer(a, j))
  names(lnapoch) <- Msumunique

  # Logarithms of pochhammer(b_i, m_i) for m_1 = 0...k, ...,  m_n = 0...k
  lnbpoch <- as.data.frame(
    matrix(nrow = max(M) + 1, ncol = n, dimnames = list(Munique, 1:n)))
  for (i in Munique) for (j in 1:n) {
    # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
    lnbpoch[i+1, j] <- lnpochhammer(b[j], i)
  }
  # matlnbpoch <- as.data.frame(matrix(nrow = max(Munique)+1, ncol = n))
  # for (i in Munique) for (j in 1:n) {
  #   # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
  #   matlnbpoch[i+1, j] <- lnpochhammer(b[j], Munique[i+1])
  # }
  # lnbpoch <- data.frame(stack(matlnbpoch), M = 0:maxM)
  # sumlnpochb <- function(ind) {
  #   # sum(mapply(function(i, j) lnbpoch[i, j], ind+1, 1:n))
  #   sum(access(ind, lnbpoch))
  # }
  # Table of all combinations of the log(pochhammer(b_i, m_i))
  gridlnbpoch <- expand.grid(as.data.frame(lnbpoch))
  # Sum of the logarithms
  sumlnpochb <- rowSums(gridlnbpoch)
  
  # Logarithms of pochhammer(c, m_1+...+m_n) for m_1 = 0...k, ...,  m_n = 0...k
  # lngpoch <- sapply(Msum, function(j) lnpochhammer(g, j))
  lngpoch <- sapply(Msumunique, function(j) lnpochhammer(g, j))
  names(lngpoch) <- Msumunique
  
  # res1 <- lnapoch + lnbpoch - lngpoch - lnfact
  # res2 <- sum(xfact * exp(res1))
  res1 <- lnapoch[Msum+1] + sumlnpochb - lngpoch[Msum+1] - sumlnfact
  # res1 <- lnapoch[Msum+1] + sumbpoch - lngpoch[Msum+1] -
  #   apply(M, 1, sumlnfact)
  res2 <- sum( prodxfact * exp(res1) )
  # res2 <- sum( apply(M, 1, prodxfact) * exp(res1) )
  
  kstep <- 5
  
  k1 <- 1:k
  result <- 0
  
  # prodxfact <- function(ind) {
  #   # prod(mapply(function(i, j) xfact[i, j], ind+1, 1:n))
  #   prod(access(ind, xfact))
  # }
  # sumlnfact <- function(ind) {
  #   # sum(mapply(function(i, j) lnfact[i, j], ind+1, 1:n))
  #   sum(access(ind, lnfact))
  # }
  # sumlnpochb <- function(ind) {
  #   # sum(mapply(function(i, j) lnbpoch[i, j], ind+1, 1:n))
  #   sum(access(ind, lnbpoch))
  # }
  
  while (abs(res2) > eps/10 & !is.nan(res2)) {
    
    epsret <- signif(abs(res2), 1)*10
    
    k <- k1[length(k1)]
    
    k1 <- k + (1:kstep)
    result <- result + res2
    
    # # M: data.frame of the indices for the nested sums
    # M <- expand.grid(rep(list(k1), n))
    # if (n > 1) {
    #   for (i in 1:(n-1)) {
    #     Mlist <- c( rep(list(0:k), n-i), rep(list(k1), i) )
    #     M <- rbind( M, expand.grid(Mlist) )
    #     for (j in 1:(n-1)) {
    #       Mlist <- Mlist[c(n, 1:(n-1))]
    #       M <- rbind(M, expand.grid(Mlist))
    #     }
    #   }
    # }
    # M1 <- M
    
    # M: data.frame of the indices for the nested sums
    M <- expand.grid(rep(list(k1), n))
    if (n > 1) {
      for (i in 1:(n-1)) {
        indsupp <- combn(n, i)
        for (j in 1:ncol(indsupp)) {
          jsupp <- indsupp[, j]
          Mlist <- vector("list", n)
          for (l in jsupp) Mlist[[l]] <- k1
          for (l in (1:n)[-jsupp]) Mlist[[l]] <- 0:k
          M <- rbind(M, expand.grid(Mlist))
        }
      }
    }
    # M2 <- M
    
    # Sum of the indices
    Msum <- rowSums(M)

    Munique <- (max(Munique)+1):max(M)
    Msumunique <- (max(Msumunique)+1):max(Msum)
    
    # Product x^{m_1} * ... * x^{m_n} for m_1, ...,  m_n given by the rows of M
    # xfact <- apply(M, 1, function(Mi) prod( x^Mi ))
    # lnfact <- apply(M, 1, function(Mi) sum(sapply(Mi, lnfactorial)))
    xfactsupp <- as.data.frame(
      matrix(nrow = length(Munique), ncol = n, dimnames = list(Munique, 1:n)))
    for (i in 1:length(Munique)) for (j in 1:n) {
      # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
      xfactsupp[i, j] <- x[j]^Munique[i]
    }
    gridxfact <- expand.grid(xfactsupp)
    for (i in 1:(n-1)) {
      indsupp <- combn(n, i)
      for (j in 1:ncol(indsupp)) {
        jsupp <- indsupp[, j]
        xfactlist <- vector("list", n)
        names(xfactlist) <- names(gridxfact)
        xfactlist[jsupp] <- xfactsupp[jsupp]
        xfactlist[-jsupp] <- xfact[-jsupp]
        gridxfact <- rbind(gridxfact, expand.grid(xfactlist))
      }
    }
    prodxfact <- apply(gridxfact, 1, prod)
    
    # Logarithm of the product m_1! * ... * m_n! for m_1, ...,  m_n given by the rows of M
    # i.e. \sum_{i=0}^n{\log{m_i!}}
    lnfactsupp <- as.data.frame(
      matrix(lfactorial(Munique), nrow = length(Munique),
             ncol = n, dimnames = list(Munique, 1:n)))
    gridlnfact <- expand.grid(lnfactsupp)
    for (i in 1:(n-1)) {
      indsupp <- combn(n, i)
      for (j in 1:ncol(indsupp)) {
        jsupp <- indsupp[, j]
        lnfactlist <- vector("list", n)
        names(lnfactlist) <- names(gridlnfact)
        lnfactlist[jsupp] <- lnfactsupp[jsupp]
        lnfactlist[-jsupp] <- lnfact[-jsupp]
        gridlnfact <- rbind(gridlnfact, expand.grid(lnfactlist))
      }
    }
    sumlnfact <- rowSums(gridlnfact)
    
    # Logarithms of pochhammer(a, m_1+...+m_n) for m_1, ...,  m_n given by the rows of M
    # lnapoch <- sapply(Msum, function(j) lnpochhammer(a, j))
    lnapochsupp <- sapply(Msumunique, function(j) lnpochhammer(a, j))
    names(lnapochsupp) <- Msumunique
    lnapoch <- c(lnapoch, lnapochsupp)
    
    # lnbpoch <- numeric(nrow(M))
    # for (i in 1:nrow(M)) {
    #   # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
    #   lnbpoch[i] <- sum( sapply(1:n, function(j) lnpochhammer(b[j], M[i, j])) )
    # }
    # Logarithms of pochhammer(b_i, m_i) for m_1, ...,  m_n given by the rows of M
    lnbpochsupp <- as.data.frame(
      matrix(nrow = length(Munique), ncol = n, dimnames = list(Munique, 1:n)))
    for (i in 1:length(Munique)) for (j in 1:n) {
      # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
      lnbpochsupp[i, j] <- lnpochhammer(b[j], Munique[i])
    }
    # Table of all combinations of the log(pochhammer(b_i, m_i))
    gridlnbpoch <- expand.grid(lnbpochsupp)
    for (i in 1:(n-1)) {
      indsupp <- combn(n, i)
      for (j in 1:ncol(indsupp)) {
        jsupp <- indsupp[, j]
        lnbpochlist <- vector("list", n)
        names(lnbpochlist) <- names(gridlnbpoch)
        lnbpochlist[jsupp] <- lnbpochsupp[jsupp]
        lnbpochlist[-jsupp] <- lnbpoch[-jsupp]
        gridlnbpoch <- rbind(gridlnbpoch, expand.grid(lnbpochlist))
      }
    }
    # Sum of the logarithms
    sumlnpochb <- rowSums(gridlnbpoch)
    
    # Logarithms of pochhammer(g, m_1+...+m_n) for m_1, ...,  m_n given by the rows of M
    # lngpoch <- sapply(Msum, function(j) lnpochhammer(g, j))
    lngpochsupp <- sapply(Msumunique, function(j) lnpochhammer(g, j))
    names(lngpochsupp) <- Msumunique
    lngpoch <- c(lngpoch, lngpochsupp)
    
    # res1 <- lnapoch + lnbpoch - lngpoch - lnfact
    # res2 <- sum(xfact * exp(res1))
    # res1 <- lnapoch[Msum+1] + apply(M, 1, sumlnpochb) - lngpoch[Msum+1] - apply(M, 1, sumlnfact)
    # res2 <- sum( apply(M, 1, prodxfact) * exp(res1) )
    res1 <- lnapoch[Msum+1] + sumlnpochb - lngpoch[Msum+1] - sumlnfact
    res2 <- sum( prodxfact * exp(res1) )
    
    # Add the new calculated values to the tables of former values
    xfact <- rbind(xfact, xfactsupp)
    lnfact <- rbind(lnfact, lnfactsupp)
    lnbpoch <- rbind(lnbpoch, lnbpochsupp)
  }
  
  result <- Re(result)
  attr(result, "epsilon") <- eps
  attr(result, "k") <- k
  
  # Returns the result of the nested sums
  return(result)
}
