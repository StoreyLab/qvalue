hist.qvalue <- function(x, ...) {
  # Histogram plot for q-value object.
  #
  # Args:
  #   x: A q-value object.
  #
  # Returns:
  #   Histogram of p-values along with lines of local FDR and q-values
  #   plotted versus p-values
  
  # Initilizations
  suppressMessages(require(ggplot2))
  suppressMessages(require(gridExtra))
  grid.newpage() 
  dfQval <- data.frame(pval = x$pvalues[order(x$pvalues)], 
                        qval = x$qvalues[order(x$pvalues)])
  dfLfdr <- data.frame(pval = x$pvalues[order(x$pvalues)],
                         qval = x$lfdr[order(x$pvalues)])
  dfQval$Variables <- "q-value"
  dfLfdr$Variables <- "local FDR"
  pi00 <- round(x$pi0, 3)
  dfc <- rbind(dfQval, dfLfdr)
  
  # Histogram figure
  p1 <-  ggplot(data.frame(pval = x$pvalues), aes(x = pval)) +
                ggtitle("p-value density histogram") +
                geom_histogram(aes(y = ..density..), colour = "black", 
                               fill = "white", binwidth = 0.04) +
                coord_cartesian(xlim = c(0, 1)) + 
                geom_abline(intercept = x$pi0, slope = 0, lty = 2, colour = "black",
                            size = 1) + 
                geom_line(data = dfc, aes(x = pval, y = qval, color = Variables)) + 
                annotate("text", label = paste("hat(pi)[0] ==", pi00), x = 0.90, 
                         y = pi00+.1, parse = TRUE, size = 3.0, colour = "black" ) +
               xlab("p-value") + 
               ylab("density")           
  print(p1, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
}
