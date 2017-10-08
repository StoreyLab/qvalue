#' @title
#' Estimate the q-values for a given set of p-values
#'
#' @description
#' Estimate the q-values for a given set of p-values.  The q-value of a
#' test measures the proportion of false positives incurred (called the
#' false discovery rate) when that particular test is called significant.
#'
#' @details
#' The function \code{\link{pi0est}} is called internally and calculates the estimate of \eqn{\pi_0}{pi_0},
#' the proportion of true null hypotheses. The function \code{\link{lfdr}} is also called internally and
#' calculates the estimated local FDR values.  Arguments for these functions can be included via \code{...} and
#' will be utilized in the internal calls made in \code{\link{qvalue}}. See \url{http://genomine.org/papers/Storey_FDR_2011.pdf}
#' for a brief introduction to FDRs and q-values.
#'
#' @param p A vector of p-values (only necessary input).
#' @param fdr.level A level at which to control the FDR. Must be in (0,1]. Optional; if this is
#' selected, a vector of TRUE and FALSE is returned that specifies
#' whether each q-value is less than fdr.level or not.
#' @param pfdr An indicator of whether it is desired to make the
#' estimate more robust for small p-values and a direct finite sample estimate of pFDR -- optional.
#' @param lfdr.out If TRUE then local false discovery rates are returned. Default is TRUE.
#' @param pi0 It is recommended to not input an estimate of pi0. Experienced users can use their own methodology to estimate
#' the proportion of true nulls or set it equal to 1 for the BH procedure.
#' @param \ldots Additional arguments passed to \code{\link{pi0est}} and \code{\link{lfdr}}.
#' 
#'
#' @return
#' A list of object type "qvalue" containing:
#' \item{call}{Function call.}
#' \item{pi0}{An estimate of the proportion of null p-values.}
#' \item{qvalues}{A vector of the estimated q-values (the main quantity of interest).}
#' \item{pvalues}{A vector of the original p-values.}
#' \item{lfdr}{A vector of the estimated local FDR values.}
#' \item{significant}{If fdr.level is specified, and indicator of whether the
#'                   q-value fell below fdr.level (taking all such q-values to be significant
#'                                                 controls FDR at level fdr.level).}
#' \item{pi0.lambda}{An estimate of the proportion of null p-values at each \eqn{\lambda}{lambda} value (see vignette).}
#' \item{lambda}{A vector of the \eqn{\lambda}{lambda} values utilized to obtain \code{pi0.lambda}.}
#'
#' @references
#' Storey JD. (2002) A direct approach to false discovery rates. Journal
#' of the Royal Statistical Society, Series B, 64: 479-498. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00346/abstract}

#' Storey JD and Tibshirani R. (2003) Statistical significance for
#' genome-wide experiments. Proceedings of the National Academy of Sciences,
#' 100: 9440-9445. \cr
#' \url{http://www.pnas.org/content/100/16/9440.full}
#'
#' Storey JD. (2003) The positive false discovery rate: A Bayesian
#' interpretation and the q-value. Annals of Statistics, 31: 2013-2035. \cr
#' \url{http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdf_1&handle=euclid.aos/1074290335}
#'
#' Storey JD, Taylor JE, and Siegmund D. (2004) Strong control,
#' conservative point estimation, and simultaneous conservative
#' consistency of false discovery rates: A unified approach. Journal of
#' the Royal Statistical Society, Series B, 66: 187-205. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2004.00439.x/abstract}
#'
#' Storey JD. (2011) False discovery rates. In \emph{International Encyclopedia of Statistical Science}. \cr
#' \url{http://genomine.org/papers/Storey_FDR_2011.pdf} \cr
#' \url{http://www.springer.com/statistics/book/978-3-642-04897-5}
#'
#' @examples
#' # import data
#' data(hedenfalk)
#' p <- hedenfalk$p
#'
#' # get q-value object
#' qobj <- qvalue(p)
#' plot(qobj)
#' hist(qobj)
#'
#' # options available
#' qobj <- qvalue(p, lambda=0.5, pfdr=TRUE)
#' qobj <- qvalue(p, fdr.level=0.05, pi0.method="bootstrap", adj=1.2)
#'
#' @author John D. Storey
#' @seealso \code{\link{pi0est}}, \code{\link{lfdr}}, \code{\link{summary.qvalue}},
#' \code{\link{plot.qvalue}}, \code{\link{hist.qvalue}}, \code{\link{write.qvalue}}
#' @keywords qvalue
#' @aliases qvalue
#' @import splines ggplot2 reshape2
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @export
qvalue <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...) {
  # Argument checks
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  }

  # Calculate pi0 estimate
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1)  {
    pi0s = list()
    pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }

  # Calculate q-value estimates
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 - (1 - p[o]) ^ m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  }
  qvals_out[rm_na] <- qvals
  # Calculate local FDR estimates
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }

  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda,
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}
