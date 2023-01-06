#' Generate m-out-of-n subsamples
#' @param m size of subsample
#' @param n size of full sample
#' @param B number of subsamples 
#' @export
#' @return A (B x m) matrix where each row is a subsample.
#' @examples
#' get.m.out.n.tuples(10, 100, 500)
get.m.out.n.tuples <- function(m, n, B) {
  n.groups <- floor(n / m)
  n.perms <- ceiling(B / n.groups)
  tuple.mat <- replicate(n.perms, {
  idx <- sample(n, n.groups * m)
      matrix(idx, ncol=m)}, simplify = FALSE)
  tuple.mat <- do.call(rbind, tuple.mat)
  tuple.mat <- tuple.mat[1:B, ]
  return(tuple.mat)
}
