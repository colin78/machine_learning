library(MCMCpack)
source("data.R")

# Improper Uniform Prior
posterior = MCMClogit(y ~ x1 + x2, data = train, seed=0)
summary(posterior)
pdf("figures/mcmc_uniform.pdf",width=7,height=5)
plot(posterior)
dev.off()
head(posterior)
hist(posterior[,1])
cov(posterior)

# Normal Prior with mean b0 = 0, precision matrix 0.001*Id
posterior2 = MCMClogit(y ~ x1 + x2, data = train, b0=0, B0=0.001, seed=0)
summary(posterior2)
pdf("figures/mcmc_normal.pdf",width=7,height=5)
plot(posterior2)
dev.off()
