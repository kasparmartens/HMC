source("densities.R")
source("HMC.R")

# create mixture of 2 gaussians
lambda = list(0.5, 0.5)
mu = list(c(-1, -1), c(2, 2))
Sigma = list(rbind(c(1, -0.5), c(-0.5, 1)), rbind(c(1, 0.75), c(0.75, 1)))


# plot density
xval = seq(-4, 4, 0.05)
yval = seq(-4, 4, 0.05)
plot(NA, xlim=range(xval), ylim=range(yval))
values = outer(xval, yval)
for(i in 1:length(xval)){
  for(j in 1:length(yval)){
    values[i, j] = gaussian_mixture_log_density(c(xval[i], yval[j]), lambda, mu, Sigma)
  }
}
contour(xval, yval, values)


# define U = - log density and gradient of U
U = function(q) (-1)*gaussian_mixture_log_density(q, lambda, mu, Sigma)
grad_U = function(q) (-1)*gaussian_mixture_log_density_and_gradient(q, lambda, mu, Sigma)$gradient

# carry out HMC sampling
res = HMC(200, U, grad_U, epsilon = 0.05, L = 10, q0 = c(-2, -2))
points(res$Q, col="red", pch=16)
