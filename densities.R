library(mvtnorm)

# evaluate mixture of MV gaussian densities at point x
gaussian_mixture_log_density = function(x, lambda, mu, Sigma){
  n_components = length(lambda)
  log_contribution = rep(NA, n_components)
  for(k in 1:n_components){
    log_contribution[k] = log(lambda[[k]]) + dmvnorm(x, mu[[k]], Sigma[[k]], log=TRUE)
  }
  log_density = log(sum(exp(log_contribution)))
  return(log_density)
}

# evaluate mixture of MV gaussian densities at point x, 
# and additionally calculate the gradient of log density
gaussian_mixture_log_density_and_gradient = function(x, lambda, mu, Sigma){
  n_components = length(lambda)
  d = length(mu[[1]])
  
  # density evaluation
  densities = rep(NA, n_components)
  for(k in 1:n_components){
    densities[k] = dmvnorm(x, mu[[k]], Sigma[[k]])
  }
  density = sum(unlist(lambda) * densities)
  
  # gradient evaluation
  grad_contribution = matrix(NA, d, n_components)
  for(k in 1:n_components){
    grad_temp = - solve(Sigma[[k]]) %*% (x - mu[[k]])
    grad_contribution[, k] = densities[k] * lambda[[k]] * grad_temp
  }
  gradient = apply(grad_contribution, 1, sum) / density
  
  return(list(log_density = log_density, gradient = gradient))
}
