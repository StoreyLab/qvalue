#' @title Histogram of p-values
#'
#' @description
#' Histogram of p-values
#'
#' @param x A q-value object.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' This function allows one to view a histogram of the p-values along with
#' line plots of the q-values and local FDR values versus p-values. The \eqn{\pi_0}{pi_0}
#' estimate is also displayed.
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
#" Storey JD, Taylor JE, and Siegmund D. (2004) Strong control,
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
#' # make histogram
#' qobj <- qvalue(p)
#' hist(qobj)
#'
#' @aliases hist, hist.qvalue
#' @author Andrew J. Bass
#' @seealso \code{\link{qvalue}}, \code{\link{plot.qvalue}}, \code{\link{summary.qvalue}}
#' @keywords histogram
#' @export
hist.qvalue <- function(x, ...) {
  pi00 <- round(x$pi0, 3)
  rm_na <- !is.na(x$pvalues)
  pvalues <- x$pvalues[rm_na]
  qvalues <- x$qvalues[rm_na]
  lfdr <- x$lfdr[rm_na]
  # Initilizations
  d <- data.frame(pvals = pvalues,
                  qvals = qvalues,
                  lfdr = lfdr,
                  pi0 = pi00)
  dm <- melt(d, id.vars = "pvals")
  # Histogram figure
  ggplot(dm, aes_string(x = 'pvals')) +
         ggtitle("p-value density histogram") +
         geom_histogram(aes_string(y = '..density..'), colour = "black",
                        fill = "white", binwidth = 0.04, center=0.02) +
         coord_cartesian(xlim = c(0, 1)) +
         geom_line(aes_string(x = 'pvals', y = 'value', color = 'variable', linetype = 'variable'), size = 1.1) +
         scale_linetype_manual(name  = "Variables",
                               values = c("lfdr"=1,"qvals"=1, "pi0"=5),
                               labels=c("q-values", "local FDR",
                                        bquote(hat(pi)[0]==.(pi00)))) +
   #   scale_size_manual(values=c(1.2, 1.2, 2.2))+
         annotate("text", label = paste("hat(pi)[0] ==", pi00), x = 0.90,
                  y = pi00 + .1, parse = TRUE, size = 3.0, colour = "black" ) +
    #     scale_colour_discrete(name  = "Variables",
    #                      breaks=c("qvals", "lfdr", "pi0"),
    #                      labels=c("q-values", "local FDR",
    #                               bquote(hat(pi)[0]==.(pi00)))) +
         xlab("p-value") +
         ylab("density")  + theme_bw() +
  scale_colour_manual(name = "Variables", values=c("red","blue","black") , breaks=c("qvals", "lfdr", "pi0"), labels=c("q-values", "local FDR",bquote(hat(pi)[0]==.(pi00))))
}
