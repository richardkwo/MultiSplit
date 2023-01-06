# Simulation for Sec 4.1.2: Dip hunting test for mixture of unit balls
# Rscript simu-unimodal-balls.R <job=1-275>

library(plyr)
library(MultiSplit)
library(future)
plan(multisession)       # run in parallel

args <- commandArgs(trailingOnly = TRUE)

n.realizations <- 20
configs.df <- expand.grid(n.splits=50,
                          n=1000,
                          B=500,
                          p=c(5, 50, 100, 500, 1000),
                          sep.size = c(0, 0.5, seq(1, 2, length.out=9)),
                          id=1:5)


if (length(args)==0) {
  job <- 101
} else {
  job <- as.integer(args[1])
}
config <- configs.df[job,]
print(config)
cat("\n")

id <- config$id
n.splits <- config$n.splits
n <- config$n
B <- config$B
p <- config$p
sep.size <- config$sep.size

simu.df <- rdply(n.realizations, {
  X <- simu.two.balls(n, p, sep=sep.size * 2 / sqrt(2 + p))
  pval.dip <- test.multisplit(X, test.unimodal.dip.hunt.single, 
                              q.null=qunif, 
                              S=list(function(x) mean(x, na.rm=TRUE), 
                                     function(x) min(x, na.rm=TRUE)), 
                              reject.larger=FALSE, 
                              n.splits=n.splits,
                              B=B,
                              verbose=TRUE)$p.value
  pval.sigclust <- sigclust::sigclust(X, nsim=1000, icovest=2)@pval
  data.frame(method=c("dip.hunting", "sigclust"), 
             pval=c(pval.dip, pval.sigclust))
}, .id = NULL)

simu.df$n.realizations <- n.realizations
simu.df$id <- id
simu.df$n.splits <- n.splits
simu.df$n <- n
simu.df$B <- B
simu.df$p <- p
simu.df$sep.size <- sep.size
simu.df$setting <- "balls"

save.filename <- tempfile("unimodal-", tmpdir="../simu/", fileext=".csv")
write.csv(simu.df,
          file=save.filename,
          row.names = FALSE)
cat(sprintf("\n * saved to %s", save.filename))