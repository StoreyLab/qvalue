#' @title Proportion of true null p-values
#' @description
#' Estimates the proportion of true null p-values, i.e., those following the Uniform(0,1) distribution.
#'
#' @param p A vector of p-values (only necessary input).
#' @param lambda The value of the tuning parameter to estimate
#' \eqn{\pi_0}{pi_0}. Must be in [0,1). Optional, see Storey (2002).
#' @param pi0.method Either "smoother" or "bootstrap"; the method for
#' automatically choosing tuning parameter in the estimation of \eqn{\pi_0}{pi_0}, 
#' the proportion of true null hypotheses.
#' @param smooth.df Number of degrees-of-freedom to use when estimating \eqn{\pi_0}{pi_0} 
#' with a smoother. Optional.
#' @param smooth.log.pi0 If TRUE and \code{pi0.method} = "smoother", \eqn{\pi_0}{pi_0} will be 
#' estimated by applying a smoother to a scatterplot of \eqn{\log(\pi_0)}{log(pi_0)} estimates 
#' against the tuning parameter \eqn{\lambda}{lambda}. Optional.
#' @param \dots Arguments passed from \code{\link{qvalue}} function.
#'
#' @details
#' If no options are selected, then the method used to estimate \eqn{\pi_0}{pi_0} is
#' the smoother method described in Storey and Tibshirani (2003). The
#' bootstrap method is described in Storey, Taylor & Siegmund (2004). A closed form solution of the
#' bootstrap method is used in the package and is significantly faster.
#'
#' @return Returns a list: 
#' \item{pi0}{A numeric that is the estimated proportion
#' of true null p-values.} 
#' \item{pi0.lambda}{A vector of the proportion of null
#' values at the \eqn{\lambda}{lambda} values (see vignette).} 
#' \item{lambda}{A vector of \eqn{\lambda}{lambda} value(s) utilized in calculating \code{pi0.lambda}.} 
#' \item{pi0.smooth}{A vector of fitted values from the
#' smoother fit to the \eqn{\pi_0}{pi_0} estimates at each \code{lambda} value
#' (pi0.method="bootstrap" returns NULL).}
#'
#' @references
#' Storey JD. (2002) A direct approach to false discovery rates. Journal
#' of the Royal Statistical Society, Series B, 64: 479-498. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00346/abstract}
#'
#' Storey JD and Tibshirani R. (2003) Statistical significance for
#' genome-wide experiments. Proceedings of the National Academy of Sciences, 
#' 100: 9440-9445. \cr
#" \url{http://www.pnas.org/content/100/16/9440.full}
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
#' # proportion of null p-values
#' nullRatio <- pi0est(p)
#' nullRatioS <- pi0est(p, lambda=seq(0.40, 0.95, 0.05), smooth.log.pi0="TRUE")
#' nullRatioM <- pi0est(p, pi0.method="bootstrap")
#' 
#' # check behavior of estimate over lambda
#' # also, pi0est arguments can be passed to qvalue
#' qobj = qvalue(p, lambda=seq(0.05, 0.95, 0.1), smooth.log.pi0="TRUE")
#' hist(qobj)
#' plot(qobj)
#' 
#' @author John D. Storey 
#' @seealso \code{\link{qvalue}}
#' @keywords pi0est, proportion true nulls
#' @aliases pi0est
#' @export
pi0est <- function(p = NULL, lambda = seq(0.05,0.95,0.05), pi0.method = c("smoother", "bootstrap"), 
                   smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
  # Check input arguments
  pi0.method = match.arg(pi0.method)
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

