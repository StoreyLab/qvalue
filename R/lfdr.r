lfdr <- function(p, pi0 = NULL, trunc = TRUE, monotone = TRUE,
                 transf = c("probit", "logit"), adj = 1.5, eps = 10 ^ -8, ...) {
  # Determines local FDR for p-values
  #
  # Args:
  #   p: A vector of p-values.
  #   pi0: Proportion of null values. Default is NULL.
  #   trunc: If TRUE, local FDR values >1 set to 1. Default is TRUE.
  #   monotone: If TRUE, preserves the order.
  #   transf: Either "probit" or "logit" transformation.
  #   adj: Numeric value that is multiple of bandwidth of distribution. Default
  #        is 1.0.
  #   eps: Numeric value that is threshold for the tails of distribution.
  #        Default is 10^-8
  #
  # Returns:
  #   A vector of local FDR values
  # Check inputs
  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range [0,1].")
  } else if (is.null(pi0)) {
    pi0 <- pi0est(p, ...)$pi0
  }
  n <- length(p)
  transf <- match.arg(transf)
  # Local FDR method for both probit and logit transformations
  if (transf == "probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)
    x <- qnorm(p)
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0 * dnorm(x) / y
  } else {
    x <- log((p + eps) / (1 - p + eps))
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1 + exp(x)) ^ 2
    lfdr <- (pi0 * dx) / y
  }
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {	
    lfdr <- lfdr[order(p)]
    for (i in 2:n) {
      if (lfdr[i] < lfdr[i - 1]) {
        lfdr[i] <- lfdr[i - 1]
      }
    }
    lfdr <- lfdr[rank(p)]
  }  
  return(lfdr)
}
