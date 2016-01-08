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
  
  # log-density evaluation
  log_contribution = rep(NA, n_components)
  for(k in 1:n_components){
    log_contribution[k] = log(lambda[[k]]) + dmvnorm(x, mu[[k]], Sigma[[k]], log=TRUE)
  }
  log_density = log(sum(exp(log_contribution)))
  
  # gradient evaluation
  grad_contribution = matrix(NA, n_components, n_components)
  for(k in 1:n_components){
    grad_temp = solve(Sigma[[k]]) %*% (x - mu[[k]])
    grad_contribution[, k] = exp(log_contribution[k]) * lambda[[k]] * grad_temp
  }
  grad_contribution = grad_contribution / log_density
  gradient = apply(grad_contribution, 1, sum)
  
  return(list(log_density = log_density, gradient = gradient))
}
