#' @title Calculate p-values from a set of observed test statistics and 
#'   simulated null test statistics
#'   
#' @description Calculates p-values from a set of observed test statistics and 
#'   simulated null test statistics
#'   
#' @param stat A vector of calculated test statistics.
#' @param stat0 A vector or matrix of simulated or data-resampled null test 
#'   statistics.
#' @param pool If FALSE, stat0 must be a matrix with the number of rows equal to
#'   the length of \code{stat}. Default is TRUE.
#'   
#' @details The argument \code{stat} must be such that the larger the value is 
#'   the more deviated (i.e., "more extreme") from the null hypothesis it is. 
#'   Examples include an F-statistic or the absolute value of a t-statistic. The
#'   argument \code{stat0} should be calculated analogously on data that 
#'   represents observations from the null hypothesis distribution. The p-values
#'   are calculated as the proportion of values from \code{stat0} that are 
#'   greater than or equal to that from \code{stat}. If \code{pool=TRUE} is 
#'   selected, then all of \code{stat0} is used in calculating the p-value for a
#'   given entry of \code{stat}. If \code{pool=FALSE}, then it is assumed that 
#'   \code{stat0} is a matrix, where \code{stat0[i,]} is used to calculate the 
#'   p-value for \code{stat[i]}. The function \code{empPvals} calculates 
#'   "pooled" p-values faster than using a for-loop.
#'   
#'   See page 18 of the Supporting Information in Storey et al. (2005) PNAS 
#'   (\url{http://www.pnas.org/content/suppl/2005/08/26/0504609102.DC1/04609SuppAppendix.pdf})
#'    for an explanation as to why calculating p-values from pooled empirical 
#'   null statistics and then estimating FDR on these p-values is equivalent to 
#'   directly thresholding the test statistics themselves and utilizing an 
#'   analogous FDR estimator.
#'   
#' @return A vector of p-values calculated as described above.
#'   
#' @references Storey JD and Tibshirani R. (2003) Statistical significance for 
#'   genome-wide experiments. Proceedings of the National Academy of Sciences, 
#'   100: 9440-9445.  \cr \url{http://www.pnas.org/content/100/16/9440.full}
#'   
#'   Storey JD, Xiao W, Leek JT, Tompkins RG, Davis RW. (2005) Significance 
#'   analysis of time course microarray experiments.  Proceedings of the
#'   National Academy of Sciences, 102 (36), 12837-12842. \cr 
#'   \url{http://www.pnas.org/content/102/36/12837.full.pdf?with-ds=yes}
#'   
#' @examples
#' # import data
#' data(hedenfalk)
#' stat <- hedenfalk$stat
#' stat0 <- hedenfalk$stat0 #vector from null distribution
#'  
#' # calculate p-values 
#' p.pooled <- empPvals(stat=stat, stat0=stat0)
#' p.testspecific <- empPvals(stat=stat, stat0=stat0, pool=FALSE)
#'  
#' # compare pooled to test-specific p-values
#' qqplot(p.pooled, p.testspecific); abline(0,1)
#' 
#' @author John D. Storey 
#' @seealso \code{\link{qvalue}}
#' @aliases empPvals
#' @keywords pvalues
#' @export
empPvals <- function(stat, stat0, pool = TRUE) {
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
    stat0 <- (stat0 - stat) >= 0
    p <- rowMeans(stat0)
    p <- pmax(p, 1 / ncol(stat0))
  }
  return(p)
}
