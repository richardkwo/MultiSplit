# dip test ---------

# for solving beta from beta = gamma^{-1}(d) ----
solve.beta <- function(d) {
  if ((6.259 <= d) & (d <= 6.29)) {
    # about 2 x pi
    return(Inf)
  } else if (d > 2 * pi) {
    gamma.fun <- function(b) {
      2 * b * beta(b-1/2, 1/2)^2
    }
    return(stats::uniroot(function(b) gamma.fun(b) - d, c(0.5, 1e3))$root)
  } else {
    gamma.fun <- function(b) {
      2^(4 * b - 1) * (b - 1) * beta(b, b)^2
    }
    return(stats::uniroot(function(b) gamma.fun(b) - d, c(1, 200))$root)
  }
}

# get reference sample -----
get.ref.dist <- function(beta.est, n, n.monte.carlo=1000) {
  if (is.infinite(beta.est)) {
    # normal
    .sample <- function() {
      stats::rnorm(n)
    }} else if (beta.est >=1) {
      # beta
      .sample <- function() {
        stats::rbeta(n, beta.est, beta.est)
      }} else {
        # rescaled student's t
        .sample <-function() {
          .df <- 2 * beta.est - 1
          stats::rt(n, df=.df, ncp=0) / sqrt(.df)
        }
      }
  replicate(n.monte.carlo, {
    z <- .sample()
    diptest::dip(z, min.is.0 = TRUE) * 2
  })
}

# estimate d -------
estimate.d <- function(x, h.sel=kedd::h.ucv) {
  h.f <- h.sel(x, deriv.order = 0)$h
  f.hat <- suppressWarnings(kedd::dkde(x, deriv.order = 0, h=h.f))
  eval.points <- f.hat$eval.points
  h.f.2nd <- h.sel(x, deriv.order=2)$h
  f.hat.2nd.deriv <- suppressWarnings(kedd::dkde(x, y=eval.points, deriv.order = 2, h=h.f.2nd))
  x0.idx <- which.max(f.hat$est.fx)
  d.hat <- abs(f.hat.2nd.deriv$est.fx[x0.idx]) / f.hat$est.fx[x0.idx]^3
  return(d.hat)
}

#' Dip test for univariate unimodalty
#'
#' @param x univariate data vector
#' @param h.sel bandwidth selector. Default to "maximum likelihood cross 
#'    validation" selector \code{\link[kedd]{h.mlcv}}.
#' @param n.monte.carlo bootstrap sample size
#' @export
#' @return p-value
#' @references 
#' Hartigan, J. A., & Hartigan, P. M. (1985). 
#' The dip test of unimodality. 
#' \emph{The Annals of Statistics}, 70-84.
#' 
#' M-Y Cheng and Peter Hall. 
#' Calibrating the excess mass and dip tests of modality. 
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 60(3):579–589, 1998.
#' @seealso This function computes a calibrated p-value. This is 
#'    different from the (conservative) p-value from \code{\link[diptest]{dip.test}}.
#' @examples
#' \dontrun{
#' dip.test.pval(rnorm(200))
#' dip.test.pval(c(rnorm(200), rnorm(200) + 3))}
dip.test.pval <- function(x, h.sel=kedd::h.mlcv, n.monte.carlo=1000) {
  stopifnot(is.vector(x))
  # normalize, somehow important for mlcv
  x <- (x - mean(x)) / stats::sd(x)
  x <- (x - min(x)) / (max(x) - min(x))
  # estimate d
  d.hat <- estimate.d(x, h.sel=h.sel)
  beta.hat <- solve.beta(d.hat)
  ref.vals <- get.ref.dist(beta.hat, length(x), n.monte.carlo = n.monte.carlo)
  dip.x <- diptest::dip(x, min.is.0 = TRUE) * 2
  pval.x <- (sum(ref.vals > dip.x) + stats::runif(1)) / (n.monte.carlo + 1)
  return(pval.x)
}

# hunt and test -------
# find direction thru k-means ------
find.direction <- function(data) {
  clustering <- flexclust::kcca(data, 2, control=list(initcent="kmeanspp"))
  dir.vec <- clustering@centers[2,] - clustering@centers[1,]
  dir.vec <- dir.vec / sqrt(sum(dir.vec^2))
  return(dir.vec)
}

#' Dip hunting test (single data split) for multivariate linear unimodality
#'
#' Test the null hypothesis that a random vector \eqn{X} comes from a unimodal 
#' distribution, specifically, \emph{linear unimodal} in the sense that 
#' every linear projection of \eqn{X} is unimodal. 
#' 
#' The test splits data into two part A and B of equal size. It runs 2-means on 
#' part A to determine a direction. Then it uses part B to conduct \code{\link{dip.test.pval}}.
#' 
#' @param X Matrix whose rows are iid copies of \eqn{X}
#' @param h.sel Bandwidth selector
#' @param n.monte.carlo bootstrap sample size
#' @return p-value
#' @export
#' @note This is a randomized test.
#' @references 
#' Hartigan, J. A., & Hartigan, P. M. (1985). 
#' The dip test of unimodality. 
#' \emph{The Annals of Statistics}, 70-84.
#' 
#' M-Y Cheng and Peter Hall. 
#' Calibrating the excess mass and dip tests of modality. 
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 60(3):579–589, 1998.
#' @seealso \code{\link{dip.test.pval}}, \code{\link{test.multisplit}}
#' @examples
#' \dontrun{
#' # unimodal
#' X <- simu.two.balls(n=1000, p=50, sep=0)
#' replicate(10, test.unimodal.dip.hunt.single(X))
#' # bimodal
#' X <- simu.two.balls(n=1000, p=50, sep=2)
#' replicate(10, test.unimodal.dip.hunt.single(X))}
test.unimodal.dip.hunt.single <- function(X, h.sel=kedd::h.mlcv, 
                                          n.monte.carlo=1000) {
  # split data
  n <- nrow(X)
  idx.1 <- sample(n, round(n/2))  # for computing the dip
  X.0 <- X[-idx.1, ]
  X.1 <- X[idx.1, ]
  # find direction
  dir.vec <- find.direction(X.0)
  # compute stat and get p-val
  y <- c(X.1 %*% dir.vec)
  dip.test.pval(y, h.sel=h.sel, n.monte.carlo = n.monte.carlo)
}


#' Simulate a mixture of two unit balls in p dimensions 
#' 
#' Draw from a 0.5-0.5 mixture of uniform distribution over unit balls that 
#' are separated at a certain distance away. 
#' 
#' @param n sample size
#' @param p dimension
#' @param sep Euclidean distance between two ball centers.
#' @return `n x p` matrix
#' @export
simu.two.balls <- function(n=100, p=3, sep=1) {
  n1 <- round(n * 0.5)
  x1 <- replicate(n1, {
    z <- stats::rnorm(p)
    z <- z / sqrt(sum(z^2))
    z * stats::rbeta(1, p, 1)
  })
  x2.mean <- stats::rnorm(p)
  x2.mean <- x2.mean / sqrt(sum(x2.mean^2)) * sep
  x2 <- replicate(n-n1, {
    z <- stats::rnorm(p)
    z <- z / sqrt(sum(z^2))
    z * stats::rbeta(1, p, 1) + x2.mean
  })
  X <- t(cbind(x1, x2))
  X[sample(nrow(X)), ]
}