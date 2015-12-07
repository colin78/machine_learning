library(MCMCpack)

source("data.R")

mcmc = MCMClogit(y ~ x1 + x2, data = train, seed=0)
summary(mcmc)
