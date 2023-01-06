#' Goodness-of-fit test of a linear quantile regression model with data splitting
#' 
#' Test the null hypothesis that `y`'s \eqn{\tau}-th quantile equals \eqn{\beta_0 + \beta^{T} X}
#' for some \eqn{\beta_0} and \eqn{\beta}. That is, this is a specification test 
#' for `quantreg::rq(y ~ ., tau=tau, data)`; see \code{\link[quantreg]{rq}}.
#' 
#' @param data data frame, whose column `y` should be the outcome
#' @param tau quantile
#' @param bsmethod bootstrap standard error estimator; see \code{\link[quantreg]{boot.rq}}
#' @export
#' @return Z-statistic that is N(0,1) distributed when model is well-specified. 
#'    When model is misspecified, the statistic tends to be a large positive value. 
#' @note This is a randomized, hunt-and-test statistic.
#' @seealso  \code{\link[quantreg]{rq}}, \code{\link{test.multisplit}}
#' @examples
#' X <- matrix(rnorm(1000 * 10), ncol=10)
#' beta <- rep(c(-1,0,1), length.out=10)
#' # well-specified
#' y <- 1 + X %*% beta + X[,1] * rt(1000, 4)
#' data <- data.frame(y=y, X=X)
#' replicate(10, gof.quantreg.test.single(data))
#' # misspecified
#' y <- 1 + X %*% beta + X[,2]^2 + X[,1] * rt(1000, 4)
#' data <- data.frame(y=y, X=X)
#' replicate(10, gof.quantreg.test.single(data))
gof.quantreg.test.single <- function(data, tau=0.5, bsmethod="xy") {
  idx.1 <- sample(nrow(data), nrow(data) / 2)
  data.1 <- data[idx.1,]
  data.2 <- data[-idx.1, ]
  # fit resid 
  fit.1 <- quantreg::rq(y ~ ., tau=tau, data=data.1, method="fn")
  data.1$y <- (stats::resid(fit.1) > 0)
  RF.pred <- ranger::ranger(y ~ ., data.1)
  # orthogonalize and fit into data.2
  data.2$r.hat <- stats::predict(RF.pred, data.2[,-1])$predictions
  ortho.fit <- stats::lm(r.hat ~ ., data=data.2[,-1])
  data.2$r.hat <- stats::resid(ortho.fit)
  # fit data.2
  fit.2 <- quantreg::rq(y ~ ., tau=tau, data=data.2, method="fn")
  se.summary <- summary(fit.2, se="boot", bsmethod=bsmethod)
  t.stat <- se.summary$coefficients["r.hat", 3]
  if (is.nan(t.stat)) {
    t.stat <- NA
  }
  return(t.stat)
}

