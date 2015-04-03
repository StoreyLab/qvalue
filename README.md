qvalue: Q-value estimation for false discovery rate control
======


The qvalue package performs false discovery rate (FDR) estimation from a collection of p-values or from a collection of test-statistics with corresponding empirical null statistics. This package produces estimates of three key quantities: q-values, the proportion of true null hypotheses (denoted by pi\_0), and local false discovery rates.

When carrying out multiple hypothesis tests, one typically starts either with a set of p-values or test-statistics.  Either quantity yields a natural ordering of tests from most significant to least significant.  For example, using p-values one would order the tests from smallest p-value (most significant) to largest p-value (least significant).  As another example, using F-statistics one would order the tests from largest F-statistic (most significant) to smallest F-statistic (least significant).

One may then ask: "If I draw a significance threshold somewhere along this list, how much can I trust the top of the list, i.e., those I choose to call statistically significant?"  Another possible question is: "Where should I draw a line of significance along this list so that we can expect that at most 10\% of the list I call significant is composed of false positives?"  We may also wish to know the reliability of a set of tests called significant for all possible thresholds simultaneously or we may want to estimate the probability that any given test is a true null hypothesis.  

The qvalue package forms various estimates that allow one to answer these and other questions.  The quantity of interest is the false discovery rate -- sometimes abbreviated as FDR -- which is roughly defined to be the expected proportion of false discoveries (also known as false positives) among all tests that are called significant.  

An overview of the FDR and its well-established methods and theory may be found in Storey (2011) (preprint freely available [here](http://genomine.org/papers/Storey_FDR_2011.pdf)).  We recommend reading the package vignette for users of qvalue who want a quick start and are unfamiliar with FDR, q-value, and local FDR estimation.


Installation and documentation
----------------------------------

To install, open R and type:

```R 
install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")
```

The vignette can be viewed by typing:

```R
browseVignettes(package = "qvalue")
```

Package overview
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
obs_stats <- hedenfalk$stat
null_stats <- hedenfalk$stat0
pvalues <- empPvals(stat = obs_stats, stat0 = null_stats)
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

For additional details, please see the package vignette.

### Point-and-click implementation
A [Shiny](http://shiny.rstudio.com "Shiny") implementation of the package written by Andrew Bass can be found [here](http://qvalue.princeton.edu "qvalue").

### Frequently asked questions
1. This package produces "adjusted p-values", so how is it possible that my adjusted p-values are smaller than my original p-values?

The q-value is not an adjusted p-value, but rather a population quantity with an explicit definition.  The package produces estimates of q-values and the local FDR, both of which are very different from p-values.  The package does not perform a Bonferroni correction on p-values, which returns "adjusted p-values" that are larger than the original p-values.  The maximum possible q-value is pi\_0, the proportion of true null hypotheses.  The maximum possible p-value is 1.  When considering a large number of hypothesis tests where there is a nontrivial fraction of true alternative p-values, we will have both an estimate pi\_0 < 1 and we will have some large p-values close to 1.   Therefore, the maximal estimated q-value will be less than or equal to the estimated pi\_0 but there will also be a number of p-values larger than the estimated pi\_0.  It must be the case then that at some point p-values become larger than estimated q-values.
