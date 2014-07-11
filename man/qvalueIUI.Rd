\name{qvalueIUI}
\alias{qvalueIUI}
\title{Interactive user interace for q-value analysis}
\description{
  Interactive feature for a q-value analysis. This function walks the user 
  through the input required for the \code{qvalue} function.
}
\usage{
qvalueIUI(outname="qvalue.obj", return.object=TRUE, ...)
}
\arguments{
  \item{outname}{Name of output object.}
  \item{return.object}{If TRUE, q-value object will be returned.}
  \item{\dots}{Arguments from \code{qvalue} function.}
}
\details{
This function walks the user through the input required for the \code{\link{qvalue}} function.  It allows the user to call \code{\link{qvalue}} using the default arguments, view/save plots, and save the results by answering a sequence of questions.  This function is automatically called when no arguments are given to the \code{\link{qvalue}} function (i.e., \code{qvalue()}).  This is intended for users with limited experience with R. 
}
\value{
  If return.object=TRUE, a qvalue object (see \code{\link{qvalue}}).
}
\references{
  Storey JD. (2002) A direct approach to false discovery rates. Journal
  of the Royal Statistical Society, Series B, 64: 479-498. \cr
  \url{http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00346/abstract}
  
  Storey JD and Tibshirani R. (2003) Statistical significance for
  genome-wide experiments. Proceedings of the National Academy of Sciences, 
  100: 9440-9445. \cr
  \url{http://www.pnas.org/content/100/16/9440.full}

  Storey JD. (2003) The positive false discovery rate: A Bayesian
  interpretation and the q-value. Annals of Statistics, 31: 2013-2035. \cr
  \url{http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdf_1&handle=euclid.aos/1074290335}

  Storey JD, Taylor JE, and Siegmund D. (2004) Strong control,
  conservative point estimation, and simultaneous conservative
  consistency of false discovery rates: A unified approach. Journal of
  the Royal Statistical Society, Series B, 66: 187-205. \cr
  \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2004.00439.x/abstract}

  Storey JD. (2011) False discovery rates. In \emph{International Encyclopedia of Statistical Science}. \cr
  \url{http://genomine.org/papers/Storey_FDR_2011.pdf} \cr
  \url{http://www.springer.com/statistics/book/978-3-642-04897-5} 
}
\author{John D. Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{qvalue}}}
\examples{
\dontrun{
#start program directly
qvalueIUI()

#starts qvalueIUI() because no arguments given to qvalue
qvalue()
}
}
\keyword{qvalue}
