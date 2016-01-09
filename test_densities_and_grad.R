source("densities.R")

# mixture of 2 gaussians
lambda = list(0.5, 0.5)
mu = list(c(-1, -1), c(2, 2))
Sigma = list(rbind(c(1, -0.5), c(-0.5, 1)), rbind(c(1, 0.75), c(0.75, 1)))

# evaluate log-density on a grid
plot(NA, xlim=c(-3, 3), ylim=c(-3, 3))
xval = seq(-3, 3, 0.05)
yval = seq(-3, 3, 0.05)
log_density = outer(xval, yval)
for(i in 1:length(xval)){
  for(j in 1:length(yval)){
    log_density[i, j] = gaussian_mixture_log_density(c(xval[i], yval[j]), lambda, mu, Sigma)
  }
}
contour(xval, yval, log_density)

# compute gradients on another grid (not so dense) and plot the arrows
x_grad = seq(-3, 3, 0.5)
y_grad = seq(-3, 3, 0.5)
for(i in 1:length(x_grad)){
  for(j in 1:length(y_grad)){
    x = x_grad[i]
    y = y_grad[j]
    obj = gaussian_mixture_log_density_and_gradient(c(x, y), lambda, mu, Sigma)
    arrows(x, y, x+obj$gradient[1]*30, y+obj$gradient[2]*30, length=0.05, col="red")
  }
}
