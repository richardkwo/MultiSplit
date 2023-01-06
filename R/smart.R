# smart aggregation -----
get.smart.agg.pval <- function(obs.vec, ref.mat, reject.larger=TRUE) {
  # obs.vec: observed agg. stat of length W
  # ref.mat: subsampling agg. stat of size B x W
  stopifnot(length(obs.vec) == ncol(ref.mat))
  W <- length(obs.vec)
  if (!reject.larger) {
    obs.vec <- -obs.vec
    ref.mat <- -ref.mat
  }
  # columnwise calibration
  ecdf.funs <- lapply(1:W, function(w) stats::ecdf(ref.mat[,w]))
  ref.mat.calibrated <- sapply(1:W, function(w) ecdf.funs[[w]](ref.mat[,w]))
  R.ref <- apply(ref.mat.calibrated, 1, max)
  R.obs <- max(sapply(1:W, function(w) ecdf.funs[[w]](obs.vec[w])))
  return(mean(R.ref > R.obs))
}