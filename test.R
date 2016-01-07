load("alphabet.RData")

# example letter
letters[["A"]]

# Draw ellipses for a mixture of MV normal distributions
draw_letter = function(symbol = "A"){
  mixture = letters[[symbol]]
  
  plot(NA, xlim = c(0, 30), ylim = c(0, 30), xlab="", ylab="")
  # choose mixture component k
  for(k in 1:mixture$n_components){
    # which isoline to draw
    for(conf_level in c(0.1, 0.5, 0.9)){
      draw_ellipse(mu = mixture$mu[[k]], Sigma = mixture$Sigma[[k]], conf_level)
    }
  }
}

# plot all the letters, one by one
symbols = names(letters)
for(i in 1:length(symbols)){
  draw_letter(symbols[i])
}
