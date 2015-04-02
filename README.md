qvalue: Q-value estimation for false discovery rate control.
======

This package takes a list of p-values resulting from the simultaneous testing of many hypotheses and estimates their q-values and local FDR values. The q-value of a test measures the proportion of false positives incurred (called the false discovery rate) when that particular test is called significant. The local FDR measures the posterior probability the null hypothesis is true given the test's p-value. Various plots are automatically generated, allowing one to make sensible significance cut-offs. Several mathematical results have recently been shown on the conservative accuracy of the estimated q-values from this software. The software can be applied to problems in genomics, brain imaging, astrophysics, and data mining.

Installation and Documentation
----------------------------------

To install, open R and type:

```R 
install.packages('devtools')
library('devtools')
install_github('jdstorey/qvalue')
```

The vignette can be viewed by typing:

```R
browseVignettes(package = "qvalue")
```

Overview of package
--------

### Functions
* `qvalue`:  Estimate the q-values for a given set of p-values.  The q-value of a test measures the proportion of false positives incurred (called the false discovery rate) when that particular test is called significant.
* `lfdr`: Estimate the local FDR values from p-values. 
* `pi0est`: Estimates the proportion of true null p-values, i.e., those following the Uniform(0,1) distribution.
* `empPvals`: Calculates p-values from a set of observed test statistics and simulated null test statistics.
* `summary`: Display summary information for a q-value object.
* `plot`: Plot of the q-value object
* `hist`: Histogram plot of the q-value object
* `write`: Write the results of the q-value object to a file.


### Quick start guide
Given a set of p-values, the qvalue object can be calculated by using the `qvalue` function:

```R
library(qvalue)
data(hedenfalk)
pvalues <- hedenfalk$p
qobj <- qvalue(p = pvalues)
```

Additionally, the qvalue object can be calculated given a set of empirical null statistics:

```R
library(qvalue)
data(hedenfalk)
obs.stat <- hedenfalk$stat
null.stat <- hedenfalk$stat0
pvalues <- empPvals(stat = obs.stat, stat0 = null.stat)
qobj <- qvalue(p = pvalues)
```

Once the qvalue object is created, estimates of the q-values, the proportion of true null hypotheses, and the local false discovery rates can be accessed from `qobj`:

```R
qvalues <- qobj$qvalues
pi0 <- qobj$pi0
lfdr <- qobj$lfdr
```

The object can be summarized and visualized by:
```R
summary(qobj)
hist(qobj)
plot(qobj)
```

For more details on how to use the package please see the package vignette.

### Point-and-click implementation
A [Shiny](http://shiny.rstudio.com "Shiny") implementation of the package written by Andrew Bass can be found [here](http://qvalue.princeton.edu "qvalue").

### Frequently asked questions
1. This package produces ``adjusted p-values'', so how is it possible that my adjusted p-values are smaller than my original p-values?
  - First, the q-value is not an "adjusted p-values", _per se_, but rather a population quantity with an explicit definition.  The package produces estimates of q-values and the local FDR, both of which are very different from p-values.  The package does not perform a Bonferroni correction on p-values, which returns "adjusted p-values" that are larger than the original p-values.  Second, the maximum possible q-value is pi\_0, the proportion of true null hypotheses.  The maximum possible p-value is 1.  When considering a large number of hypothesis tests where there is a nontrivial fraction of true alternative p-values, we will have both an estimate pi\_0 < 1 and we will have some large p-values close to 1.   Therefore, the maximal estimated q-value will be less than or equal to the estimated pi\_0 but there will also be a number of p-values larger than the estimated pi\_0.  It must be the case then that at some point p-values become larger than estimated q-values.  This is expected and it is not a bug.
