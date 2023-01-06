# Simulation for Sec 4.1.2: Dip hunting test for mixture of multivariate t's
# Rscript simu-unimodal-t.R <job=1-275>

library(plyr)
library(MASS)
library(MultiSplit)
library(future)
plan(multisession)       # run in parallel

args <- commandArgs(trailingOnly = TRUE)

n.realizations <- 20
configs.df <- expand.grid(n.splits=50,
                          n=1000,
                          B=5,
                          rho=0.5, 
                          nu=4,
                          p=c(5, 50, 100, 500, 1000),
                          sep.size = c(0, seq(0.5, 2, length.out=11)),
                          id=1:5)

if (length(args)==0) {
  job <- 101
} else {
  job <- as.integer(args[1])
}

simu.mv.t <- function(n=100, p=3, rho=0.5, nu=4) {
  Sigma <- rho^abs(outer(1:p, 1:p, FUN="-"))
  X <- mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma)
  u <- rgamma(n, nu/2, nu/2)
  return(diag(1/u) %*% X)
}


simu.two.mv.ts <- function(n=100, p=3, rho=0.5, nu=4, sep.size=1) {
  n1 <- round(n * 0.5)
  Sigma <- rho^abs(outer(1:p, 1:p, FUN="-"))
  mu <- eigen(Sigma, symmetric = TRUE)$vectors[,2] * sqrt(p) * sep.size
  X.1 <- simu.mv.t(n=n1, p=p, rho=rho, nu=nu)
  X.2 <- simu.mv.t(n=n-n1, p=p, rho=rho, nu=nu)
  X.2 <- X.2 + matrix(rep(mu, nrow(X.2)), nrow=nrow(X.2), byrow = TRUE)
  X <- rbind(X.1, X.2)
  X <- X[sample(nrow(X)), ]
  return(X)
}


config <- configs.df[job,]
print(config)
cat("\n")

id <- config$id
n.splits <- config$n.splits
n <- config$n
B <- config$B
p <- config$p
nu <- config$nu
rho <- config$rho
sep.size <- config$sep.size

simu.df <- rdply(n.realizations, {
  X <- simu.two.mv.ts(n=n, p=p, sep.size=sep.size, rho=rho, nu=nu)
  pval.dip <- test.multisplit(X, test.unimodal.dip.hunt.single, 
                              q.null=qunif, 
                              S=list(function(x) mean(x, na.rm=TRUE), 
                                     function(x) min(x, na.rm=TRUE)), 
                              reject.larger=FALSE, 
                              n.splits=n.splits,
                              B=B,
                              verbose=TRUE)$p.value
  pval.sigclust <- sigclust::sigclust(X, nsim=1000, icovest=2)@pval
  data.frame(method=c(names(pval.dip), "sigclust"),
             pval=c(pval.dip, pval.sigclust))
}, .id = "realization")

simu.df$id <- id
simu.df$n.splits <- n.splits
simu.df$n <- n
simu.df$B <- B
simu.df$p <- p
simu.df$nu <- nu
simu.df$rho <- rho
simu.df$sep.size <- sep.size
simu.df$setting <- "t"


save.filename <- tempfile("unimodal-", tmpdir="../simu/", fileext=".csv")
write.csv(simu.df,
          file=save.filename,
          row.names = FALSE)
cat(sprintf("\n * saved to %s", save.filename))