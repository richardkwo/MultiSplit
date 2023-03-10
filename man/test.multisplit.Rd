% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank.R
\name{test.multisplit}
\alias{test.multisplit}
\title{Test with multiple splits}
\usage{
test.multisplit(
  data,
  test.single,
  q.null = stats::qnorm,
  S = function(x) mean(x, na.rm = TRUE),
  reject.larger = TRUE,
  n.splits = 50,
  B = NULL,
  packages = c(),
  verbose = FALSE
)
}
\arguments{
\item{data}{a data frame or matrix with iid rows}

\item{test.single}{a function that computes a randomized (e.g., through data
splitting) test statistic from \code{data}. Calling \code{test.single(data)}
should return a real number and should handle \code{data} with a different
number of rows.}

\item{q.null}{Quantile function of the asymtotic null distribution of \code{test.single},
e.g., \code{qnorm} for a Z-statistic or \code{qunif} for a p-value. \strong{Note}:
quantile functions other than \code{qnorm} and \code{qunif} may not be theoretically
supported.}

\item{S}{aggregation function(s). It can be either a single aggregation function
(e.g., \code{S=mean}) or a list of aggregation functions (e.g., \code{S=c(mean, max)}).
For the latter case, a p-value that automatically adapts to the most
powerful aggregation function in the list will be returned.}

\item{reject.larger}{The side of the test. Use \code{FALSE} if \code{test.single} returns
a p-value.}

\item{n.splits}{Number of applying \code{test.single} to \code{data} (e.g., number of
data splits)}

\item{B}{Number of subsamples (similar to number of bootstrap replications).
When set to \code{NULL}, will set to default value \code{80 * log(nrow(data))}.}

\item{packages}{A list of packages that \code{test.single} depends on to be passed
to \code{\link[future]{future}}. Default to \code{c()}.}

\item{verbose}{If \code{TRUE}, will print details.}
}
\value{
A list consists of the aggregated observed statistic(s) and its
p-value.
}
\description{
Compute a randomized test statistic multiple times (e.g., through multiple
data splits), aggregate them with function \code{S} and compute the p-value
of the aggregated statistic.
}
\note{
Package \pkg{future} is used to compute subsampling in parallel. To
enable parallelization, the user is responsible to select the right "plan"
with \code{\link[future]{plan}} beforehand.
}
\examples{
\dontrun{
# test multivariate unimodality
library(future)
plan(multisession)   # run in parallel
# unimodal
X <-simu.two.balls(n=1000, p=50, sep=0)
test.multisplit(X, test.unimodal.dip.hunt.single, 
                q.null=qunif, 
                S=list(function(x) mean(x, na.rm=TRUE), 
                function(x) min(x, na.rm=TRUE)), 
                reject.larger=FALSE, 
                n.splits=10,
                verbose=TRUE)
# bimodal
X <-simu.two.balls(n=1000, p=50, sep=1)
test.multisplit(X, test.unimodal.dip.hunt.single, 
                q.null=qunif, 
                S=list(function(x) mean(x, na.rm=TRUE), 
                function(x) min(x, na.rm=TRUE)), 
                reject.larger=FALSE, 
                n.splits=10,
                verbose=TRUE)                
}
}
