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


textsize = 18


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

mu = w_iter[2:3,]
sigma <- as.matrix(V[2:3,2:3]) # matrix(c(1, 0.6, 0.6, 1), nrow=2)
rng   <- 15
N     <- 50
X     <- seq(mu[1]-rng*sigma[1, 1], mu[1]+rng*sigma[1, 1], length.out=N)
Y     <- seq(mu[2]-rng*sigma[2, 2], mu[2]+rng*sigma[2, 2], length.out=N)
set.seed(123)
library(mvtnorm)
genZ <- function(x, y) { dmvnorm(cbind(x, y), mu, sigma) }
matZ <- outer(X, Y, FUN="genZ")
c = contour(X, Y, matZ, nlevels=5, col="blue", drawlabels=FALSE,
            main="Contours for 2D-normal density",
            xlab=expression(x[1]), ylab=expression(x[2]), lwd=2)
pdf("figures/mcmc_uniform_2d_2.pdf", width=8, height=5)
c
dev.off()

library(mvtnorm)
library(ggplot2)
set.seed(1)
n <- 1e3
c95 <- qchisq(.95, df=2)
c90 = qchisq(.9, df=2)
c65 = qchisq(.65, df=2)
rho <- 0.8  #correlation 
Sigma <- matrix(c(1, rho, rho, 1), 2, 2) # Covariance matrix
x <- rmvnorm(n, mean=c(0, 0), Sigma)
z  <- p95 <- rep(NA, n)
for(i in 1:n){
  z[i] <- x[i, ] %*% solve(Sigma, x[i, ])
  p95[i] <- (z[i] < c95)
}
data <- data.frame(x, z, p95)
p <- ggplot(data, aes(X1, X2)) + geom_point(aes(colour = p95))
p + stat_ellipse(type = "norm")
p <- ggplot(data, aes(X1, X2)) + geom_point(aes(colour = p95))
p + stat_ellipse(type = "norm")

ellips <- function(center = c(0,0), c=c95, rho=-0.8, npoints = 100){
  t <- seq(0, 2*pi, len=npoints)
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  a <- sqrt(c*eigen(Sigma)$values[2])
  b <- sqrt(c*eigen(Sigma)$values[1])
  x <- center[1] + a*cos(t)
  y <- center[2] + b*sin(t)
  X <- cbind(x, y)
  R <- eigen(Sigma)$vectors
  data.frame(X%*%R)
}
dat <- ellips(center=c(0, 0), c=c95, rho, npoints=100)
dat90 <- ellips(center=c(0, 0), c=c90, rho, npoints=100)
dat65 <- ellips(center=c(0, 0), c=c65, rho, npoints=100)

p + geom_path(data=dat, aes(x=X1, y=X2), colour='blue') + 
 geom_path(data=dat90, aes(x=X1, y=X2), colour='blue') + 
 geom_path(data=dat65, aes(x=X1, y=X2), colour='blue')


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
n = 1000
# Simulate 10000 points from multivariate normal N(w, V) and N(w_iter, V_iter)
# (excludes w_0 component)
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
