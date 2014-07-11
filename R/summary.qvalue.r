summary.qvalue <- function (object, cuts = c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1), 
                      digits = getOption("digits"), ...) {
  # Provides p value, q value and local FDR results under various significant
  # levels
  #
  # Args:
  #   qobj: A q-value object returned by the Qvalue function.
  #   cuts: Various cutoffs for p, q, local FDR values.
  #   digits: Specifies the length of values for output.
  #
  # Returns:
  #   Prints a table of counts for the the number of features under certain
  #   significant levels.
  # Call
  cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  # Proportion of nulls
  cat("pi0:", format(object$pi0, digits=digits), "\n", sep="\t")
  cat("\n")
  # Number of significant values for p-value, q-value and local FDR
  cat("Cumulative number of significant calls:\n")
  cat("\n")
  counts <- sapply(cuts, function(x) c("p-value"=sum(object$pvalues < x),
                  "q-value"=sum(object$qvalues < x), 
                  "local FDR"=sum(object$lfdr < x)))
  colnames(counts) <- paste("<", cuts, sep="")
  print(counts)
  cat("\n")
}
