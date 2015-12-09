library(MCMCpack)
source("data.R")

# Improper Uniform Prior
posterior = MCMClogit(y ~ x1 + x2, data = train, seed=0)
summary(posterior)
pdf("figures/mcmc_uniform.pdf",width=7,height=5)
plot(posterior)
dev.off()
posterior = data.frame(posterior)
tail(posterior)
hist(posterior[,1])
cov(posterior)

# Normal Prior with mean b0 = 0, precision matrix 0.001*Id
posterior2 = MCMClogit(y ~ x1 + x2, data = train, b0=0, B0=0.001, seed=0)
summary(posterior2)
pdf("figures/mcmc_normal.pdf",width=7,height=5)
plot(posterior2)
dev.off()
posterior2 = data.frame(posterior2)

# Obtain vb_logit_fit, vb_logit_fit_iter results
w = read.csv("results/w.csv", header=FALSE)
w_iter = read.csv("results/w_iter.csv", header=FALSE)
V = read.csv("results/V.csv", header=FALSE)
V_iter = read.csv("results/V_iter.csv", header=FALSE)

# Simulate 10000 points from multivariate normal N(w, V) and N(w_iter, V_iter)
# (excludes w_0 component)
textsize = 18
n = 100000
samples = mvrnorm(n, w[2:3,],V[2:3,2:3])
samples_iter = mvrnorm(n, w_iter[2:3,],V_iter[2:3,2:3])
samples = data.frame(samples)
names(samples) = c("w1", "w2")
samples_iter = data.frame(samples_iter)
names(samples_iter) = c("w1", "w2")
g2 = ggplot(data=samples, aes(x=w1, y=w2)) + geom_point() +
  theme(text=element_text(size=textsize, family="CM Sans")) + 
  xlab(expression(w[1])) + ylab(expression(w[2]))
g2 + geom_density2d()

samples[,1]

# Plot MCMC sample results
library(ggplot2)
# # see https://github.com/wch/fontcm
# install.packages("extrafont")
library(extrafont)
# font_install('fontcm')
loadfonts()


library(MASS)
indep = posterior[seq(1,10000,10),]
indep2 = posterior2[seq(1,10000,10),]


g = ggplot(data=indep, aes(x=x1, y=x2)) + geom_point() +
  theme(text=element_text(size=textsize, family="CM Sans")) + 
  xlab(expression(w[1])) + ylab(expression(w[2]))
g = g + geom_density2d()
ggsave("figures/mcmc_uniform_2d.pdf", g, width=8, height=5)

mu    <- c(1, 3)
sigma <- matrix(c(1, 0.6, 0.6, 1), nrow=2)
rng   <- 2.5
N     <- 50
X     <- seq(mu[1]-rng*sigma[1, 1], mu[1]+rng*sigma[1, 1], length.out=N)
Y     <- seq(mu[2]-rng*sigma[2, 2], mu[2]+rng*sigma[2, 2], length.out=N)
set.seed(123)
library(mvtnorm)
genZ <- function(x, y) { dmvnorm(cbind(x, y), mu, sigma) }
matZ <- outer(X, Y, FUN="genZ")
c = contour(X, Y, matZ, main="Contours for 2D-normal density")
ggsave("figures/mcmc_uniform_2d_2.pdf", c, width=8, height=5)

pdf("figures/mcmc_uniform_2d_2.pdf", width=8, height=5)
c
g
dev.off()


require(mvtnorm)
x1 = seq(-3, 3, length.out=200)
x2 = seq(-3, 3, length.out=200)
z = matrix(0, length(x1), length(x2))
for (i in 1:length(x1)) {
  a = x1
  b = x2[i]
  z[,i] = dmvnorm(cbind(a,b))
}
image(x1,x2,z)
# reshape the data
# require(reshape2)
# dat <- melt(z)
z
# use geom_raster to mimic image

gg <- ggplot(as.data.frame(z), aes(x=Var2, y=Var1)) + geom_point()
gg = gg + geom_contour(aes())
gg
gg <- gg + geom_raster()
gg <- gg + coord_equal()
gg <- gg + scale_fill_gradient(low="red", high="yellow")
gg <- gg + scale_x_continuous(expand = c(0, 0))
gg <- gg + scale_y_continuous(expand = c(0, 0))
gg <- gg + theme_bw()
gg
