% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unimodal.R
\name{simu.two.balls}
\alias{simu.two.balls}
\title{Simulate a mixture of two unit balls in p dimensions}
\usage{
simu.two.balls(n = 100, p = 3, sep = 1)
}
\arguments{
\item{n}{sample size}

\item{p}{dimension}

\item{sep}{Euclidean distance between two ball centers.}
}
\value{
\verb{n x p} matrix
}
\description{
Draw from a 0.5-0.5 mixture of uniform distribution over unit balls that
are separated at a certain distance away.
}
