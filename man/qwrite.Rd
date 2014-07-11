\name{qwrite}
\alias{qwrite}
\title{Write results to file}
\description{
  Write the results of the q-value object to a file. 
}
\usage{
qwrite(qobj, filename = "my-qvalue-results.txt", ...)
}
\arguments{
  \item{qobj}{A q-value object.}
  \item{filename}{Output filename (optional).}
  \item{\dots}{Additional arguments; currently unused.}

}
\details{
  The output file lists the estimate of \eqn{\pi_0}{pi_0}, which is the
  proportion of true null hypotheses. It also lists each p-value and
  corresponding q-values and local FDR values, one per line. If an FDR significance level was 
  specified in the call to \code{\link{qvalue}}, the significance level 
  is printed below the estimate of \eqn{\pi_0}{pi_0}, and an indicator 
  of significance is included as a fourth column for each p-value and q-value.
}
\value{
  Nothing of interest.
}
\author{John D. Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{qvalue}}, \code{\link{plot.qvalue}}, \code{\link{summary.qvalue}}}
\examples{
# import data
data(hedenfalk)
p <- hedenfalk$p

# write q-value object
qobj <- qvalue(p)
qwrite(qobj, file="myresults.txt")
}
\keyword{write}
\keyword{qwrite}
