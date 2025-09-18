# smart aggregation -----
get.smart.agg.pval <- function(obs.vec, ref.mat, side) {
  # obs.vec: observed agg. stat of length W
  # ref.mat: subsampling agg. stat of size B x W
  stopifnot(side %in% c(-1,0,1))
  stopifnot(length(obs.vec) == ncol(ref.mat))
  W <- length(obs.vec)
  if (side==-1) {
    obs.vec <- -obs.vec
    ref.mat <- -ref.mat
  } else if (side==0) {
    obs.vec <- abs(obs.vec)
    ref.mat <- abs(ref.mat)
  }
  # columnwise calibration
  ecdf.funs <- lapply(1:W, function(w) stats::ecdf(ref.mat[,w]))
  ref.mat.calibrated <- sapply(1:W, function(w) ecdf.funs[[w]](ref.mat[,w]))
  R.ref <- apply(ref.mat.calibrated, 1, max)
  R.obs <- max(sapply(1:W, function(w) ecdf.funs[[w]](obs.vec[w])))
  return(list(R.ref, R.obs))
}