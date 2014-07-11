pi0est <- function(p = NULL, lambda = seq(0.05,0.95,0.05), pi0.method = "smoother", 
                   smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
  # Estimates the proportion of null p-values, pi0
  #
  # Args:
  #   p: A vector of p-values.
  #   lambda: A vector of the tuning parameters to estimate pi0.
  #   pi0.method: Either "smoother" or "bootstrap"; the method for 
  #               automatically choosing tuning parameter in the estimation of
  #               pi0, the proportion of true null hypotheses.
  #   smooth.df: Degrees of freedom to use in smoother; Default is 3.
  #   smooth.log.pi0: If TRUE, smoothing is done on a log scale; Default is 
  #                   FALSE.
  #
  # Returns:
  #   pi0: An estimate of the proportion of null p-values.
  #   lambda: A vector of lambda value(s) choosen.
  #   pi0.lambda: A vector of the proportion of null values at lambda.
  #   pi0.smooth: A vector of pi0 estimates at each lambda for smoother
  #               method. If bootstrap method, value is NULL.
  # Check input arguments
  m <- length(p)
  lambda <- sort(lambda) # guard against user input
  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  } else if (ll > 1 && ll < 4) {
    stop(cat("ERROR:", paste("length(lambda)=", ll, ".", sep=""), 
             "If length of lambda greater than 1,", 
             "you need at least 4 values."))
  } else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  # Determines pi0
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL 
  } else {
    pi0 <- sapply(lambda, function(l) mean(p >= l) / (1 - l))
    pi0.lambda <- pi0
    # Smoother method approximation
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      } else {    
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    } else if (pi0.method == "bootstrap") {  
      # Bootstrap method closed form solution by David Robinson 
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W / (m ^ 2 * (1 - lambda) ^ 2)) * (1 - W / m) + (pi0 - minpi0) ^ 2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    } else {
      stop('ERROR: pi0.method must be one of "smoother" or "bootstrap".')
    }
  }  
  if (pi0 <= 0) {
    stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda, 
              lambda = lambda, pi0.smooth = pi0Smooth))
}

