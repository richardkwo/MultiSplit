#' Test with multiple splits
#'
#' Compute a randomized test statistic multiple times (e.g., through multiple
#' data splits), aggregate them with function \code{S} and compute the p-value
#' of the aggregated statistic.
#' 
#' @param data a data frame or matrix with iid rows
#' @param test.single a function that computes a randomized (e.g., through data 
#'    splitting) test statistic from `data`. Calling `test.single(data)` 
#'    should return a real number and should handle `data` with a different 
#'    number of rows. 
#' @param q.null Quantile function of the asymtotic null distribution of `test.single`,
#'    e.g., `qnorm` for a Z-statistic or `qunif` for a p-value. \strong{Note}: 
#'    quantile functions other than `qnorm` and `qunif` may not be theoretically
#'    supported.
#' @param S aggregation function(s). It can be either a single aggregation function
#'    (e.g., `S=mean`) or a list of aggregation functions (e.g., `S=c(mean, max)`).
#'    For the latter case, a p-value that automatically adapts to the most 
#'    powerful aggregation function in the list will be returned.
#' @param reject.larger The side of the test. Use `FALSE` if `test.single` returns
#'    a p-value.
#' @param n.splits Number of applying `test.single` to `data` (e.g., number of 
#'    data splits)
#' @param B Number of subsamples (similar to number of bootstrap replications).
#'    When set to `NULL`, will set to default value `80 * log(nrow(data))`.
#' @param packages A list of packages that `test.single` depends on to be passed 
#'    to \code{\link[future]{future}}. Default to `c()`.
#' @param verbose If `TRUE`, will print details. 
#' @export
#' @return A list consists of the aggregated observed statistic(s) and its
#'    p-value.
#' @note Package \pkg{future} is used to compute subsampling in parallel. To
#'    enable parallelization, the user is responsible to select the right "plan" 
#'    with \code{\link[future]{plan}} beforehand.
#' @examples
#' \dontrun{
#' # test multivariate unimodality
#' library(future)
#' plan(multisession)   # run in parallel
#' # unimodal
#' X <-simu.two.balls(n=1000, p=50, sep=0)
#' test.multisplit(X, test.unimodal.dip.hunt.single, 
#'                 q.null=qunif, 
#'                 S=list(function(x) mean(x, na.rm=TRUE), 
#'                 function(x) min(x, na.rm=TRUE)), 
#'                 reject.larger=FALSE, 
#'                 n.splits=10,
#'                 verbose=TRUE)
#' # bimodal
#' X <-simu.two.balls(n=1000, p=50, sep=1)
#' test.multisplit(X, test.unimodal.dip.hunt.single, 
#'                 q.null=qunif, 
#'                 S=list(function(x) mean(x, na.rm=TRUE), 
#'                 function(x) min(x, na.rm=TRUE)), 
#'                 reject.larger=FALSE, 
#'                 n.splits=10,
#'                 verbose=TRUE)                
#' }
test.multisplit <- function(data, test.single, 
                            q.null=stats::qnorm,
                            S=function(x) mean(x, na.rm=TRUE),
                            reject.larger=TRUE,
                            n.splits=50, B=NULL, 
                            packages=c(), 
                            verbose=FALSE) {
  n <- nrow(data)
  if (is.null(B)) {
    B <- floor(5 * n / log(n))
  } else {
    B <- as.integer(B)
  }
  stopifnot(is.function(q.null))
  # subsampling
  m <- round(n / log(n))
  tuple.mat <- get.m.out.n.tuples(m, n, B)
  if (verbose) {
    message(sprintf("Subsampling with %d tuples (m = %d) and %d splits ...", B, m, n.splits))
  }
  pval.jobs <- plyr::llply(1:B, function(b) {
    future::future({
      data.sub <- data[tuple.mat[b,], ]  # of size m
      # iterate over splits
      replicate(n.splits, {
        test.single(data.sub)
      })
    }, packages = packages, seed=TRUE)
  }, .progress = ifelse(verbose, "text", "none"))
  # retrieve value
  T.mat <- plyr::laply(pval.jobs, future::value)
  T.mat.transformed <- (rank(T.mat, ties.method = "random") - 1/2) / length(T.mat)  # rank transform
  T.mat.transformed <- matrix(T.mat.transformed, nrow=B)
  # apply inverse CDF
  T.mat.transformed <- q.null(T.mat.transformed)
  if (verbose) {
    message(sprintf("Max change from rank transform = %g", max(abs(T.mat.transformed - T.mat))))
  }
  # observed T
  T.obs.vec <- replicate(n.splits, test.single(data))
  # aggregation 
  if (!is.list(S)) {
    # single aggregation function ------
    if (verbose) {
      message(sprintf("Single aggregation function (reject for %s values)", ifelse(reject.larger, "larger", "smaller")))
    }
    T.obs <- S(T.obs.vec)
    T.sub <- apply(T.mat.transformed, 1, S)
    if (reject.larger) {
      p.value <- mean(T.sub > T.obs)
    } else {
      p.value <- mean(T.sub < T.obs)
    }
  } else {
    # multiple aggregation functions -----
    if (verbose) {
      message(sprintf("Adapting to the best among %d aggregation functions (reject for %s values)", length(S), ifelse(reject.larger, "larger", "smaller")))
    }
    T.obs <- sapply(S, function(.s) .s(T.obs.vec))
    T.sub <- sapply(S, function(.s) apply(T.mat.transformed, 1, .s))
    p.value <- get.smart.agg.pval(T.obs, T.sub, reject.larger=reject.larger)
  }
  return(list(p.value=p.value, T.obs=T.obs))
}
