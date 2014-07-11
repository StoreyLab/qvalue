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
  suppressMessages(require(ggplot2))
  suppressMessages(require(gridExtra))
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
    p1.smooth <- geom_line(data = spi0.df, aes(x = lambda, y = pi0), colour="red")
  }        
  # Subplots
  p1 <-  ggplot(pi0.df, aes(x = lambda, y = pi0)) +
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
                          aes(x = pvalue, y = qvalue)) + 
               xlab("p-value") +
               ylab("q-value")  + 
               geom_line()
  p3 <- ggplot(data.frame(qCuttOff = q.ord[q.ord >= rng[1] & q.ord <= rng[2]],
                          sig=(1 + sum(q.ord < rng[1])):sum(q.ord <= rng[2])), 
                          aes(x = qCuttOff, y = sig)) + 
               xlab("q-value cut-off") + 
               ylab("significant tests")  + 
               geom_line()
  p4 <- ggplot(data.frame(sig = (1 + sum(q.ord < rng[1])):sum(q.ord <= rng[2]),
                          expFP = q.ord[q.ord >= rng[1] & q.ord <= rng[2]] * 
                                (1 + sum(q.ord < rng[1])):sum(q.ord <= rng[2])),
                          aes(x = sig, y = expFP)) + 
               xlab("significant tests") + 
               ylab("expected false positives")  + 
               geom_line()
  multPlot <- arrangeGrob(p1, p2, p3, p4)
  return(print(multPlot))
}
