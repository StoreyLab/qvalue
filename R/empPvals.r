empPvals <- function(stat, stat0, pool = TRUE) {
  # Calculates p-values from a matrix of t values of the alternative and null
  # distributions
  #
  # Args:
  #   stat: A vector of calculated test statistics.
  #   stat0: A vector or matrix of simulated or data-resampled null test 
  #          statistics. 
  #   pool: Alternative way to generate p-values. Default is TRUE.
  #
  # Returns:
  #   p values for the measured distribution
  m <- length(stat)
  n <- ncol(stat0)
  # Calculates p-values
  if (pool == TRUE) {
    if (is.matrix(stat0)) {
      stat0 <- as.vector(stat0)
    }
    m0 <- length(stat0) 
    v <- c(rep(TRUE, m), rep(FALSE, m0))
    v <- v[order(c(stat, stat0), decreasing = TRUE)]
    u <- 1:length(v)
    w <- 1:m
    p <- (u[v == TRUE] - w) / m0
    p <- p[rank(-stat)]
    p <- pmax(p, 1/m0)
  } else {
    if (is.vector(stat0)) {
      stop("stat0 must be a matrix.")
    } else if (n == m) {
      stat0 <- t(stat0)
    } else if (nrow(stat0) != m){
      stop("Number of rows of stat0 must equal length of stat.")
    }
    stat0 <- (stat0 - matrix(stat, nrow = m, ncol = n)) >= 0
    p <- rowMeans(stat0)
    p <- pmax(p, 1 / ncol(stat0))
  }
  return(p)
}