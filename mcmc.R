library(MCMCpack)

source("data.R")

posterior = MCMClogit(y ~ x1 + x2, data = train, b0=0, B0=0.001, seed=0)
plot(posterior)
summary(posterior)
