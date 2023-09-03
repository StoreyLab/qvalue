#' Q-value on truncated p-values
#' 
#' Calculate \code{\link{qvalue}}, but with p-values re-scaled to be in `[0, 1]` range 
#' 
#' @seealso \code{\link{qvalue}}
#' 
#' @export
qvalue_truncp <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...) {
  p <- p / max(p)
  qvalue(p, fdr.level, pfdr, lfdr.out, pi0, ...)
}

