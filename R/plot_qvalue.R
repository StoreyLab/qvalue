#' @title Plotting function for q-value object
#' @description
#'  Graphical display of the q-value object
#'
#'  @param x A q-value object.
#'  @param rng Range of q-values to show. Optional
#'  @param \ldots Additional arguments. Currently unused.
#'
#' @details
#' The function plot allows one to view several plots:
#' \enumerate{
#'  \item The estimated \eqn{\pi_0}{pi_0} versus the tuning parameter
#'  \eqn{\lambda}{lambda}.
#'  \item The q-values versus the p-values.
#'  \item The number of significant tests versus each q-value cutoff.
#'  \item The number of expected false positives versus the number of
#'  significant tests.
#'  }
#'  
#' This function makes four plots. The first is a plot of the
#' estimate of \eqn{\pi_0}{pi_0} versus its tuning parameter
#' \eqn{\lambda}{lambda}. In most cases, as \eqn{\lambda}{lambda} 
#' gets larger, the bias of the estimate decreases, yet the variance 
#' increases. Various methods exist for balancing this bias-variance 
#' trade-off (Storey 2002, Storey & Tibshirani 2003, Storey, Taylor 
#' & Siegmund 2004). Comparing your estimate of \eqn{\pi_0}{pi_0} to this 
#' plot allows one to guage its quality. The remaining three plots
#' show how many tests are called significant and how many false
#' positives to expect for each q-value cut-off. A thorough discussion of
#' these plots can be found in Storey & Tibshirani (2003).
#'
#' @return
#' Nothing of interest.
#'
#' @references
#' Storey JD. (2002) A direct approach to false discovery rates. Journal
#' of the Royal Statistical Society, Series B, 64: 479-498. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00346/abstract}
#'
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
#" \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2004.00439.x/abstract}
#'
#' Storey JD. (2011) False discovery rates. In \emph{International Encyclopedia of Statistical Science}. \cr
#' \url{http://genomine.org/papers/Storey_FDR_2011.pdf} \cr
#' \url{http://www.springer.com/statistics/book/978-3-642-04897-5} 
#'
#' @examples
#' # import data
#' data(hedenfalk)
#' p <- hedenfalk$p
#' qobj <- qvalue(p) 
#'  
#" # view plots for q-values between 0 and 0.3:
#' plot(qobj, rng=c(0.0, 0.3))
#'
#' @author John D. Storey, Andrew J. Bass
#' @seealso \code{\link{qvalue}}, \code{\link{write.qvalue}}, \code{\link{summary.qvalue}}
#' @keywords plot
#' @aliases plot, plot.qvalue
#' @export
plot.qvalue <- function(x, rng = c(0.0, 0.1), ...) { 
  # Plotting function for q-object.
  #
  # Args:
  #   x: A q-value object returned by the qvalue function.
  #   rng: The range of q-values to be plotted (optional).
  #
  # Returns
  #   Four plots-
  #     Upper-left: pi0.hat(lambda) versus lambda
  #     Upper-right: q-values versus p-values
  #     Lower-left: number of significant tests per each q-value cut-off
  #     Lower-right: number of expected false positives versus number of
  #                  significant tests 
  # Initilizations
  plot.call <- match.call()
  q.ord <- x$qval[order(x$pval)]
  if (min(q.ord) > rng[2]) {
    rng <- c(min(q.ord), quantile(q.ord, 0.1))
  }
  p.ord <- x$pval[order(x$pval)]
  lambda <- x$lambda
  pi0Smooth <- x$pi0.smooth
  if (length(lambda) == 1) {
    lambda <- sort(unique(c(lambda, seq(0, max(0.90, lambda), 0.05))))
  }
  pi0 <- x$pi0.lambda  
  pi00 <- round(x$pi0, 3)
  pi0.df <- data.frame(lambda = lambda, pi0 = pi0)
  # Spline fit- pi0Smooth NULL implies bootstrap
  if (is.null(pi0Smooth)) {
    p1.smooth <- NULL
  } else {
    spi0.df <- data.frame(lambda = lambda, pi0 = pi0Smooth)
    p1.smooth <- geom_line(data = spi0.df, aes_string(x = 'lambda', y = 'pi0'), 
                           colour="red")
  }        
  # Subplots
  p1 <-  ggplot(pi0.df, aes_string(x = 'lambda', y = 'pi0')) +
                geom_point() + 
                p1.smooth +
                geom_abline(intercept = pi00,
                            slope = 0, 
                            lty = 2, 
                            colour = "red",
                            size = .6) + 
                xlab(expression(lambda)) + 
                ylab(expression(hat(pi)[0](lambda))) + 
                xlim(min(lambda) - .05, max(lambda) + 0.05) +
                annotate("text", label = paste("hat(pi)[0] ==", pi00),
				                 x = min(lambda, 1) + (max(lambda) - min(lambda))/20, 
				                 y = x$pi0 - (max(pi0) - min(pi0))/20, 
                         parse = TRUE, size = 3)
  p2 <- ggplot(data.frame(pvalue = p.ord[q.ord >= rng[1] & q.ord <= rng[2]],
                          qvalue = q.ord[q.ord >= rng[1] & q.ord <= rng[2]]), 
                          aes_string(x = 'pvalue', y = 'qvalue')) + 
               xlab("p-value") +
               ylab("q-value")  + 
               geom_line()
  p3 <- ggplot(data.frame(qCuttOff = q.ord[q.ord >= rng[1] & q.ord <= rng[2]],
                          sig=(1 + sum(q.ord < rng[1])):sum(q.ord <= rng[2])), 
                          aes_string(x = 'qCuttOff', y = 'sig')) + 
               xlab("q-value cut-off") + 
               ylab("significant tests")  + 
               geom_line()
  p4 <- ggplot(data.frame(sig = (1 + sum(q.ord < rng[1])):sum(q.ord <= rng[2]),
                          expFP = q.ord[q.ord >= rng[1] & q.ord <= rng[2]] * 
                                (1 + sum(q.ord < rng[1])):sum(q.ord <= rng[2])),
                          aes_string(x = 'sig', y = 'expFP')) + 
               xlab("significant tests") + 
               ylab("expected false positives")  + 
               geom_line()
  multiplot(p1, p2, p3, p4, cols = 2)
}
