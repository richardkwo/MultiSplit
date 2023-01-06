# Simulation for Sec 4.1.3: Goodness-of-fit testing of quantile regression
# Rscript simu-gof.R <job=1-360>

library(MASS)
library(quantreg)
library(plyr)
library(MultiSplit)
library(future)
plan(multisession)       # run in parallel

# configs ------
args <- commandArgs(trailingOnly = TRUE)
n.realizations <- 20
configs.df <- expand.grid(n=c(1000),
                           p=c(10, 20),
                           tau=c(0.5),
                           v=c(0, 1, 2, 3, 4, 5),
                           nonlinear.fun=c(1,2),
                           err=c("t", "exp"),
                           spec=c(1,2),
                           id=1:5)
configs.df <- configs.df[!(configs.df$err=="t" & configs.df$spec==2), ]

if (length(args)==0) {
  job <- 101
} else {
  job <- as.integer(args[1])
}

n.splits <- 50
B <- 500

config <- configs.df[job,]
print(config)
cat("\n")

id <- config$id
n <- config$n
p <- config$p
v <- config$v
tau <- config$tau
nonlinear.fun <- config$nonlinear.fun
err <- config$err
spec <- config$spec

# data generation ------
simu.data <- function(n=1000, p=7, v=0, nonlinear.fun=1, err="t", spec=1) {
  Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
  X <- mvrnorm(n, mu=rep(0, p), Sigma=Sigma)
  beta.0 <- rep_len(c(-1,2,0), length.out = p)
  # error
  stopifnot(err %in% c("t", "exp", "normal"))
  if (err=="t") {
    eps <- rt(n, df=3)
  } else if (err=="exp") {
    eps <- rexp(n)
  } else if (err=="normal") {
    eps <- rnorm(n)
  }
  # non-linearity
  stopifnot(nonlinear.fun %in% c(1,2))
  if (nonlinear.fun==1) {
    nl <- c(sqrt(rowSums(X[,1:2]^2))) * 4
  } else if (nonlinear.fun==2) {
    nl <- exp(-(1 + X[,2] + X[,3])) * 2
  }
  # generate
  stopifnot(spec %in% c(1,2))
  if (spec==1) {
    y <- 1 + X %*% beta.0 + (1 + X[,2] + X[,3]) * eps + v / sqrt(n) * nl
  } else {
    y <- 1 + X %*% beta.0 + (1 + X[,2] + X[,3] + v / sqrt(n) * nl) * eps
  }
  data <- as.data.frame(cbind(y,X))
  colnames(data) <- c("y", paste("x", 1:p, sep="."))
  return(data)
}

# simulation ------
simu.df <- rdply(n.realizations, {
  data <- simu.data(n=n, p=p, v=v, nonlinear.fun=nonlinear.fun, err=err, spec=spec)
  
  pval <- test.multisplit(data, gof.quantreg.test.single, 
                              q.null=qnorm, 
                              S=function(x) mean(x, na.rm=TRUE), 
                              reject.larger=TRUE, 
                              n.splits=n.splits,
                              B=B,
                              verbose=TRUE)$p.value
  
  data.frame(method="hunt.and.test", pval=pval)
}, .id = "realization")

simu.df$n.splits <- n.splits
simu.df$B <- B
simu.df$id <- id
simu.df$n <- n
simu.df$p <- p
simu.df$v <- v
simu.df$tau <- tau
simu.df$nonlinear.fun <- nonlinear.fun
simu.df$err <- err
simu.df$spec <- spec

save.filename <- tempfile("gof-", tmpdir="../simu/", fileext=".csv")
write.csv(simu.df,
          file=save.filename,
          row.names = FALSE)
cat(sprintf("\n * saved to %s", save.filename))

