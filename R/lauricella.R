lauricella <- function(a, b, g, x, eps = 1e-06) {
  #' Lauricella Function
  #'
  #' Computes the Lauricella function.
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
  #' @details If \eqn{n} is the length of the \eqn{b} and \code{x} vectors, the Lauricella function is given by:
  #' 
  #' \eqn{F_D^{(n)}\left(a; b_1, ..., b_n; g; x_1, ..., x_n\right) = \sum_{m_1 \geq 0} ... \sum_{m_n \geq 0}{ \frac{ (a)_{m_1+...+m_n}(b_1)_{m_1} ... (b_n)_{m_n} }{ (g)_{m_1+...+m_n} } \frac{x_1^{m_1}}{m_1!} ... \frac{x_1^{m_n}}{m_n!} } }
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
  #' @author Nizar Bouhlel, Pierre Santagostini
  #' @references N. Bouhlel, A. Dziri, Jullback-Leibler Divergence Between Multivariate Generalized Gaussian Distributions.
  #' IEEE Signal Processing Letters, vol. 26 no. 7, July 2019.
  #' \doi{10.1109/LSP.2019.2915000}
  #' @export
  
  # Number of variables
  n <- length(x)
  
  # Do x and b have the same length?
  if (length(b) != n)
    stop("x and b must have the same length")
  
  # Condition for the convergence: are all abs(x) < 1 ?
  if (any(abs(x) >= 1))
    stop("The series does not converge for these x values.")
  
  k <- 5
  
  # M: data.frame of the indices for the nested sums
  # (i.e. all arrangements of n elements from {0:k})
  M <- expand.grid(rep(list(0:k), n))
  
  Msum <- apply(M, 1, sum)
  
  # Product x^{m_1}/m_1! * ... * x^{m_n}/m_n! for m_1 = 0...k, ...,  m_n = 0...k
  xfact <- apply(M, 1, function(Mi) prod( x^Mi / factorial(Mi) ))
  
  # pochhammer(a, m_1+...+m_n) for m_1 = 0...k, ...,  m_p = 0...k
  apoch <- sapply(Msum, function(j) pochhammer(a, j))
  
  bpoch <- numeric(nrow(M))
  for (i in 1:nrow(M)) {
    # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
    bpoch[i] <- prod( sapply(1:n, function(j) pochhammer(b[j], M[i, j])) )
  }
  
  # pochhammer(c, m_1+...+m_n) for m_1 = 0...k, ...,  m_p = 0...k
  gpoch <- sapply(Msum, function(j) pochhammer(g, j))
  
  res1 <- sum( xfact * apoch * bpoch / gpoch )
  
  kstep <- 5
  
  k1 <- 1:k
  result <- 0
  
  while (abs(res1) > eps/10 & !is.nan(res1)) {
    
    epsret <- signif(abs(res1), 1)*10
    
    k <- k1[length(k1)]
    
    k1 <- k + (1:kstep)
    result <- result + res1
    
    # M: data.frame of the indices for the nested sums
    M <- as.data.frame(matrix(nrow = 0, ncol = n))
    if (n > 1) {
      for (i in 1:(n-1)) {
        Mlist <- c( rep(list(0:k), n-i), rep(list(k1), i) )
        M <- rbind( M, expand.grid(Mlist) )
        for (j in 1:(n-1)) {
          Mlist <- Mlist[c(n, 1:(n-1))]
          M <- rbind(M, expand.grid(Mlist))
        }
      }
    }
    M <- rbind( M, expand.grid(rep(list(k1), n)) )
    
    Msum <- apply(M, 1, sum)
    
    # Product x^{m_1}/m_1! * ... * x^{m_n}/m_n! for m_1, ...,  m_n given by the rows of M
    xfact <- apply(M, 1, function(Mi) prod( x^Mi / factorial(Mi) ))
    
    # pochhammer(a, m_1+...+m_n) for m_1, ...,  m_n given by the rows of M
    apoch <- sapply(Msum, function(j) pochhammer(a, j))
    
    bpoch <- numeric(nrow(M))
    for (i in 1:nrow(M)) {
      # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
      bpoch[i] <- prod( sapply(1:n, function(j) pochhammer(b[j], M[i, j])) )
    }
    
    # pochhammer(c, m_1+...+m_n) for m_1 = 0...k, ...,  m_p = 0...k
    gpoch <- sapply(Msum, function(j) pochhammer(g, j))
    
    res1 <- sum( xfact * apoch * bpoch / gpoch )
    
  }
  
  if (is.nan(res1)) {
    
    res1 <- 0
    k1 <- k
    
    while(!is.nan(res1)) {
      
      if (res1 > 0)
        epsret <- abs(res1)*10
      
      k <- k1
      
      k1 <- k + 1
      result <- result + res1
      
      # M: data.frame of the indices for the nested sums
      M <- as.data.frame(matrix(nrow = 0, ncol = n))
      if (n > 1) {
        for (i in 1:(n-1)) {
          Mlist <- c( rep(list(0:k), n-i), rep(list(k1), i) )
          M <- rbind( M, expand.grid(Mlist) )
          for (j in 1:(n-1)) {
            Mlist <- Mlist[c(n, 1:(n-1))]
            M <- rbind(M, expand.grid(Mlist))
          }
        }
      }
      M <- rbind( M, rep(k1, n) )
      
      Msum <- apply(M, 1, sum)
      
      # Product x^{m_1} * ... * x^{m_n} for m_1 = 0...k, ...,  m_n = 0...k
      xfact <- apply(M, 1, function(Mi) prod( x^Mi / factorial(Mi) ))
      
      # pochhammer(a, m_1+...+m_n) for m_1 = 0...k, ...,  m_p = 0...k
      apoch <- sapply(Msum, function(j) pochhammer(a, j))
      
      bpoch <- numeric(nrow(M))
      for (i in 1:nrow(M)) {
        # Product pochhammer(b_1,m_1) * ...* pochhammer(b_n, m_n)
        bpoch[i] <- prod( sapply(1:n, function(j) pochhammer(b[j], M[i, j])) )
      }
      
      # pochhammer(c, m_1+...+m_n) for m_1 = 0...k, ...,  m_p = 0...k
      gpoch <- sapply(Msum, function(j) pochhammer(g, j))
      
      res1 <- sum( xfact * apoch * bpoch / gpoch )
      
    }
  }
  
  if (is.nan(res1)) {
    epsret <- signif(epsret, 1)
    warning("Cannot reach the precision ", eps, " due to NaN\n",
            "Number of iteration: ", k, "\n",
            "Precision reached:", epsret)
    attr(result, "epsilon") <- epsret
  } else {
    # result <- result + res1
    attr(result, "epsilon") <- eps
  }
  
  attr(result, "k") <- k
  
  # Returns the result of the nested sums
  return(result)
}
