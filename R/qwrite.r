qwrite <- function(qobj, filename="my-qvalue-results.txt", ...) {
  # Writes results to a file
  #
  # Args:
  #   qobj: A q-value object returned by the qvalue function.
  #   filename: The name of the file where the results are written.
  #
  # Returns:
  #   A file sent to "filename" with the following:
  #   First row: The estimate of the proportion of true negatives, pi0.
  #   Second row: FDR significance level (if specified).
  #   Third row and below: The p-values (1st column), the estimated q-values 
  #                        (2nd column), the estimate local FDR values (3rd column)
  #                        and indicator of significance level
  #                        if appropriate (4th column).
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
