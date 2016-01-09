HMC = function(n_samples, U, grad_U, epsilon, L, q0){
  d = length(q0)
  Q = matrix(NA, n_samples, d)
  accepted = rep(NA, n_samples)
  Q[1, ] = q0
  for(i in 2:n_samples){
    obj = HMC_one_step(U, grad_U, epsilon, L, current_q = Q[i-1, ])
    Q[i, ] = obj$q
    accepted[i] = obj$accepted
  }
  return(list(Q = Q, acceptance_rate = mean(accepted, na.rm=T)))
}

# the following originates mainly from http://www.cs.toronto.edu/~radford/ham-mcmc-simple
HMC_one_step = function(U, grad_U, epsilon, L, current_q){
  q = current_q
  p = rnorm(length(q), 0, 1) # momentum
  current_p = p
  
  # Make a half step for momentum at the beginning
  
  p = p - epsilon * grad_U(q) / 2
  
  # Alternate full steps for position and momentum
  
  for (i in 1:L){
    # Make a full step for the position
    
    q = q + epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory
    
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  
  # Make a half step for momentum at the end.
  
  p = p - epsilon * grad_U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
    return (list(q = q, accepted = TRUE))  # accept
  }
  else{
    return (list(q = current_q, accepted = FALSE))  # reject
  }
}
