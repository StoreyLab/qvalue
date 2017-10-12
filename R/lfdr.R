#' @title Estimate local False Discovery Rate (FDR)
#'
#' @description
#' Estimate the local FDR values from p-values.
#'
#' @param p A vector of p-values (only necessary input).
#' @param pi0 Estimated proportion of true null p-values. If NULL, then \code{\link{pi0est}} is called.
#' @param trunc If TRUE, local FDR values >1 are set to 1. Default is TRUE.
#' @param monotone If TRUE, local FDR values are non-decreasing with increasing p-values. Default is TRUE; this is recommended.
#' @param transf Either a "probit" or "logit" transformation is applied to the p-values so that a local FDR estimate can be formed that
#' does not involve edge effects of the [0,1] interval in which the p-values lie.
#' @param adj Numeric value that is applied as a multiple of the smoothing bandwidth used in the density estimation. Default is \code{adj=1.0}.
#' @param eps Numeric value that is threshold for the tails of the empirical p-value distribution. Default is 10^-8.
#' @param \ldots Additional arguments, passed to \code{\link{pi0est}}.
#'
#' @details It is assumed that null p-values follow a Uniform(0,1) distribution.
#'   The estimated proportion of true null hypotheses \eqn{\hat{\pi}_0}{pi_0} is either
#'   a user-provided value or the value calculated via \code{\link{pi0est}}.
#'   This function works by forming an estimate of the marginal density of the
#'   observed p-values, say \eqn{\hat{f}(p)}{f(p)}. Then the local FDR is estimated as
#'   \eqn{{\rm lFDR}(p) = \hat{\pi}_0/\hat{f}(p)}{lFDR(p) = pi_0/f(p)}, with
#'   adjustments for monotonicity and to guarantee that \eqn{{\rm lFDR}(p) \leq
#'   1}{lFDR(p) <= 1}.  See the Storey (2011) reference below for a concise
#'   mathematical definition of local FDR.
#'
#' @return
#' A vector of estimated local FDR values, with each entry corresponding to the entries of the input p-value vector \code{p}.
#'
#' @references
#' Efron B, Tibshirani R, Storey JD, and Tisher V. (2001) Empirical Bayes analysis
#' of a microarray experiment. Journal of the American Statistical Association, 96: 1151-1160. \cr
#' \url{http://www.tandfonline.com/doi/abs/10.1198/016214501753382129}
#'
#' Storey JD. (2003) The positive false discovery rate: A Bayesian
#' interpretation and the q-value. Annals of Statistics, 31: 2013-2035. \cr
#' \url{http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdf_1&handle=euclid.aos/1074290335}
#'
#' Storey JD. (2011) False discovery rates. In \emph{International Encyclopedia of Statistical Science}. \cr
#' \url{http://genomine.org/papers/Storey_FDR_2011.pdf} \cr
#' \url{http://www.springer.com/statistics/book/978-3-642-04897-5}
#'
#' @examples
#' # import data
#' data(hedenfalk)
#' p <- hedenfalk$p
#' lfdrVals <- lfdr(p)
#'
#' # plot local FDR values
#' qobj = qvalue(p)
#' hist(qobj)
#'
#' @author John D. Storey
#' @seealso \code{\link{qvalue}}, \code{\link{pi0est}}, \code{\link{hist.qvalue}}
#' @aliases lfdr
#' @keywords local False Discovery Rate, lfdr
#' @export
lfdr <- function(p, pi0 = NULL, trunc = TRUE, monotone = TRUE,
                 transf = c("probit", "logit"), adj = 1.5, eps = 10 ^ -8, ...) {
  # Check inputs
  lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
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
    o <- order(p, decreasing = FALSE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  return(lfdr_out)
}
