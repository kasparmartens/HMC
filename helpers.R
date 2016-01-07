# helper function for drawing an ellipse
draw_ellipse = function(mu, Sigma, conf_level = 0.95, ...){
  theta = seq(0, 2*pi, 0.01)
  
  eig = eigen(Sigma)
  lambdas = eig$values
  Gamma = eig$vectors
  
  const = qchisq(conf_level, 2)
  z1 = sqrt(lambdas[1] * const) * cos(theta)
  z2 = sqrt(lambdas[2] * const) * sin(theta)
  Z = rbind(z1, z2)
  # Alternatiiv oleks:
  # Z = sqrt(const) * diag(sqrt(lambdas)) %*% rbind(cos(theta), sin(theta))
  X = mu + Gamma %*% Z
  lines(X[1, ], X[2, ], ...)
}
