qvalue <- function(p = NULL, fdr.level = NULL, pfdr = FALSE, ...) {
  # Estimating the q-values for a given set of 
  # p-values or t-values. The methodology comes from a series of papers on 
  # false discovery rates by John D. Storey.
  #
  # Args:
  #   p: A vector of p-values.
  #   fdr.level: A level at which to control the FDR.
  #   lambda: The value of the tuning parameter to estimate pi0.
  #   stat: A matrix of t values.
  #   stat0: A matrix of random t values representative of the null statistic. 
  #   pi0.method: Either "smoother" or "bootstrap"; the method for 
  #               automatically choosing tuning parameter in the estimation of
  #               pi0, the proportion of true null hypotheses.
  #   pfdr: If TRUE, an indicator of whether it is desired to make the estimate
  #         more robust for small p-values and a direct finite sample estimate 
  #         of pFDR; Default is FALSE.
  #   smooth.df: Degrees of freedom to use in smoother; Default is 3.
  #   smooth.log.pi0: If TRUE, smoothing is done on a log scale; Default is 
  #                   FALSE. 
  #
  # Returns:
  #   call: Gives the function call.
  #   pi0: An estimate of the proportion of null p-values.
  #   qvalues: A vector of the estimated q-values.
  #   pvalues: A vector of the original p-values.
  #   lfdr: A vector of the local FDR values.
  #   significant: If fdr.level is TRUE, an indicator of whether the q-values
  #                fell below fdr.level (taking all such q-values to be 
  #                significant controls FDR at level fdr.level).
  #   lambda: A vector of lambda value(s) choosen.
  #   pi0.lambda: A vector of the proportion of null values at each lambda.

  # Argument checks
  if (is.null(p)) {
    qvalueIUI()
    return(invisible(NULL))
  } else if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  } 
  
  # Calculate pi0 estimate
  pi0s <- pi0est(p, ...)
  
  # Calculate q-value estimates
  m <- length(p)
  u <- order(p)
  v <- rank(p, ties.method="max") 
  if (pfdr) {
    qvals <- (pi0s$pi0 * m * p) / (v * (1 - (1 - p) ^ m))
  } else {
    qvals <- (pi0s$pi0 * m * p) / v
  }
  qvals[u[m]] <- min(qvals[u[m]], 1)
  for (i in (m - 1):1) {
    qvals[u[i]] <- min(qvals[u[i]], qvals[u[i+1]])
  }

  # Calculate local FDR estimates  
  lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
  
  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals,
                   pvalues = p, lfdr = lfdr, fdr.level = fdr.level, 
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda, 
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals, 
                   pvalues = p, lfdr = lfdr, pi0.lambda = pi0s$pi0.lambda, 
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}
