#' @title Write results to file
#' @description
#' Write the results of the q-value object to a file. 
#' @param qobj A q-value object.
#' @param filename Output filename (optional).
#' @param \dots Additional arguments; currently unused.
#'
#' @details
#' The output file lists the estimate of \eqn{\pi_0}{pi_0}, which is the
#' proportion of true null hypotheses. It also lists each p-value and
#' corresponding q-values and local FDR values, one per line. If an FDR significance level was 
#' specified in the call to \code{\link{qvalue}}, the significance level 
#' is printed below the estimate of \eqn{\pi_0}{pi_0}, and an indicator 
#' of significance is included as a fourth column for each p-value and q-value.
#'
#' @return 
#'  Nothing of interest.
#'
#' @examples 
#' # import data
#' data(hedenfalk)
#' p <- hedenfalk$p
#'  
#' # write q-value object
#  qobj <- qvalue(p)
#' qwrite(qobj, file="myresults.txt")
#' @author John D. Storey \email{jstorey@@princeton.edu}, Andrew J. Bass
#' @seealso \code{\link{qvalue}}, \code{\link{plot.qvalue}}, \code{\link{summary.qvalue}}
#' @aliases qwrite
#' @keywords write, qwrite
#' @export
qwrite <- function(qobj, filename="my-qvalue-results.txt", ...) {
  if (any(names(qobj) == "fdr.level")) {
    cat(c("# pi0:", "\t", qobj$pi0, "\t", " ", "\t", " ", "\n#\n"), 
        file=filename,
        append=FALSE)
    cat(c("# FDR level:", "\t", qobj$fdr.level, "\t", " ", "\t", " ", "\n#\n"), 
        file=filename,
        append=TRUE)
    cat(c("p-value", "\t", "q-value", "\t", "lfdr", "\t", "significant", "\n"), 
        file=filename, 
        append=TRUE)
    write(t(cbind(qobj$pval, qobj$qval, qobj$lfdr, qobj$significant)),
          file=filename,
          ncolumns=4, 
          append=TRUE, 
          sep="\t")
  } else {
    cat(c("# pi0:", "\t", qobj$pi0, "\t", " ", "\n#\n"),
        file=filename,
        append=FALSE)
    cat(c("p-value", "\t", "q-value", "\t", "lfdr", "\n"),
        file=filename, 
        append=TRUE)
    write(t(cbind(qobj$pval, qobj$qval, qobj$lfdr)), 
          file=filename,
          ncolumns=3,
          append=TRUE, 
          sep="\t")
  }
}
