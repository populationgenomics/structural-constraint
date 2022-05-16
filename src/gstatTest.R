
# LT 8/05

# some tests on kriging and spatial covariance functions
# using package gstat

# objective:
# 1. test applicability to the structural constraint project (simulator component)
# 2. test performance for large models (n>1e6): application to CellRegMap

require(gstat)
require(MASS)

# specify a simple exponential kernel
m = vgm(1, 'Exp', 10, 0)

# sample size for exp,
N = 100

# generate some random 1D coordinates
x = rnorm(N,1,2)

# dependent variable
oe = rnorm(100,1)

obs = data.frame(x, oe)

#coordinates(obs)=~x
g <- gstat(formula = oe~1, data = obs, model = m, locations = ~x)

# cannot be done in one dimension
obs$y = rnorm(N, 0,2)

g <- gstat(formula = oe~1, data = obs, model = m, locations = ~x+y)
plot(variogram(g))

predict(g, obs, nsim=1)

v = variogram(oe~1, locations = ~x+y, obs)
fit.variogram(v, vgm(0, "Exp", NA, NA))

# flat variogram, as expected

# generate correlated data
dmat = as.matrix(dist(obs[,c('x', 'y')]))

# exponential spatial covariance with range l
sigma=exp(-dmat/0.5)

obs$oe2 = mvrnorm(1, mu=rep(0,100),Sigma=sigma)
v2 = variogram(oe2~1, locations = ~x+y, obs)
plot(v2)
fit.variogram(v2, vgm(NA, "Exp", NA))

# same in 3d
obs3d = data.frame(x = rnorm(100, sd=2), y = rnorm(100, sd=2), z = rnorm(100, sd=2))
dmat = as.matrix(dist(obs3d[,c('x', 'y', 'z')]))

# exponential spatial covariance with range l
sigma=exp(-dmat/0.5)
obs3d$oe = mvrnorm(1, mu=rep(0,100),Sigma=sigma)
v2 = variogram(oe~1, locations = ~x+y+z, obs3d)
plot(v2)
fit.variogram(v2, vgm(NA, "Exp", NA))

# fitted parameters match true parameters
# This confirms proper construction of spatial covariances

