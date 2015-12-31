#' @title Write results to file
#' @description Write the results of the q-value object to a file.
#' @param x A q-value object.
#' @param file Output filename (optional).
#' @param sep Separation between columns.
#' @param eol Character to print at the end of each line.
#' @param na String to use when there are missing values.
#' @param row.names logical. Specify whether row names are to be printed.
#' @param col.names logical. Specify whether column names are to be printed.
#'
#' @details The output file includes: (i) p-values, (ii)
#' q-values (iii) local FDR values, and (iv) the estimate of \eqn{\pi_0}{pi_0},
#' one per line. If an FDR significance
#' level was specified in the call to \code{\link{qvalue}}, the significance
#' level is printed and an indicator of significance is included.
#'
#' @return Nothing of interest.
#'
#' @examples
#' # import data
#' data(hedenfalk)
#' p <- hedenfalk$p
#'
#' # write q-value object
#' qobj <- qvalue(p)
#' write.qvalue(qobj, file="myresults.txt")
#' @author John D. Storey, Andrew J. Bass
#' @seealso \code{\link{qvalue}}, \code{\link{plot.qvalue}},
#'   \code{\link{summary.qvalue}}
#' @aliases write.qvalue
#' @keywords write.qvalue
#' @export
write.qvalue <- function(x, file = NULL, sep = " ", eol = "\n", na = "NA",
                         row.names = FALSE, col.names = TRUE) {
  if (class(x) != "qvalue") {
    stop("x must be a qvalue object.")
  }
  d <- data.frame(pvalue = x$pval,
                  qvalue = x$qval,
                  lfdr = x$lfdr,
                  pi0 = x$pi0)
  if (any(names(x) == "fdr.level")) {
    d$significant <- x$significant
    d$fdr.level <- x$fdr.level
    write.table(as.matrix(d), file = file, sep = sep, eol = eol, na = na,
                row.names = row.names, col.names = col.names)
  } else {
    write.table(as.matrix(d), file = file, sep = sep, eol = eol, na = na,
                row.names = row.names, col.names = col.names)
  }
}
