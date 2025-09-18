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
#'    (e.g., `S=mean`) or a list of aggregation functions (e.g., `S=list(mean, max)`).
#'    For the latter case, a p-value that automatically adapts to the most 
#'    powerful aggregation function in the list will be returned.
#' @param side The side of the test: `1` for right-sided (reject larger values of the 
#'    test statistic), `-1` for left-sided (reject smaller values) and `0` for two-sided. 
#'    Use `-1` if `test.single` returns a p-value.
#' @param n.splits Number of applying `test.single` to `data` (e.g., number of 
#'    data splits)
#' @param B Number of subsamples (similar to number of bootstrap replications).
#'    When set to `NULL`, will set to default value `80 * log(nrow(data))`.
#' @param packages A list of packages that `test.single` depends on to be passed 
#'    to \code{\link[future]{future}}. Default to `c()`.
#' @param verbose If `TRUE`, will print details. 
#' @param .plot If `TRUE`, will plot the observed value along with the subsamples. 
#'    Default to `FALSE`.
#' @export
#' @return A list consists of the following fields: 
#'   - `p.value`: p-value (without smoothing the subsample's empirical distribution)
#'   - `p.value.kde`: p-value (with Gaussian kernel smoothing of the subsample's empirical distribution)
#'   - `T.obs.vec`: the vector of observed statistics resulting from `n.splits` splits
#'   - `T.agg.obs`: the aggregated observed statistic (a vector if `S` contains more than one function)
#'   - `T.agg.sub`: the subsample counterparts of the aggregated observed statistic (a matrix 
#'              if `S` contains more than one function, whose columns correspond
#'              to the functions in `S`)
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
                            side=1,
                            n.splits=50, B=NULL, 
                            packages=c(), 
                            verbose=FALSE, 
                            .plot=FALSE) {
  stopifnot(side %in% c(-1,0,1))
  if (side==1) {
    .side <- "right"
  } else if (side==-1) {
    .side <- "left"
  } else {
    .side <- "two"
  }
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
      message(sprintf("Single aggregation function (%s-sided test)", .side))
    }
    T.obs <- S(T.obs.vec)
    T.sub <- apply(T.mat.transformed, 1, S)
    } else {
      # multiple aggregation functions -----
      if (verbose) {
        message(sprintf("Adapting to the best among %d aggregation functions (%s-sided test)", length(S), .side))
      }
      T.obs.multiple <- sapply(S, function(.s) .s(T.obs.vec))
      T.sub.multiple <- sapply(S, function(.s) apply(T.mat.transformed, 1, .s))
      R.list <- get.smart.agg.pval(T.obs.multiple, T.sub.multiple, side)
      T.sub <- R.list[[1]]
      T.obs <- R.list[[2]]
      # stats have been transformed to right sided
      side <- 1
    }
  # estimate the cdf and get p-value --------
  .cdf.est <- estimate.cdf(T.sub)
  .cdf.smoothed <- .cdf.est$.cdf.smoothed
  .ecdf <- .cdf.est$.ecdf
  .x <- .cdf.est$.x
  if (side==1) {
    p.value <- 1 - .ecdf(T.obs)
    p.value.kde <- 1 - .cdf.smoothed(T.obs)
  } else if (side==-1) {
    p.value <- .ecdf(T.obs)
    p.value.kde <- .cdf.smoothed(T.obs)
  } else {
    p.value <- .ecdf(-abs(T.obs)) + 1 - .ecdf(abs(T.obs))
    p.value.kde <- .cdf.smoothed(-abs(T.obs)) + 1 - .cdf.smoothed(abs(T.obs))
  }
  if (.plot) {
    graphics::split.screen(c(2,1))
    graphics::screen(1)
    plot(.x, .ecdf(.x), type="l", 
         xlab="stat", ylab="ECDF", 
         main=sprintf("p-value=%f", p.value))
    graphics::rug(T.sub, ticksize=0.08, lwd=0.7)
    graphics::abline(v=T.obs, col="red", lwd=1.5)
    graphics::screen(2)
    plot(.x, .cdf.smoothed(.x), type="l",
         xlab="stat", ylab="CDF (ks)", 
         main=sprintf("p-value=%f", p.value.kde))
    graphics::rug(T.sub, ticksize=0.08, lwd=0.7)
    graphics::abline(v=T.obs, col="red", lwd=1.5)
    graphics::close.screen(all = TRUE)
  }
  result <- list(p.value=p.value, 
                 p.value.kde=p.value.kde, 
                 T.obs.vec=T.obs.vec)
  if (!is.list(S)) {
    result$T.agg.obs <- T.obs
    result$T.agg.sub <- T.sub
  } else {
    result$T.agg.obs <- T.obs.multiple
    result$T.agg.sub <- T.sub.multiple
  }
  return(result)
}
