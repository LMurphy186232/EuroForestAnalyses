###############################################################################
# Disperse shell
# 
# Filling out the pieces needed for the disperse model...
# 
# Lora Murphy, 5/20/2020
###############################################################################

#-----------------------------------------------------------------------#
# Lognormal isotropic model
#-----------------------------------------------------------------------#
lognormal.isotropic <- function(TLP, alpha, X0, Xb) {
  incr <- 0.25
  maxdist <- 250
  
  #########################
  # Calculate normalizer
  #########################
  # Distance increments
  dist <- seq(from=incr, to=maxdist, by=incr)
  
  # Midpoints - halfway between the increments
  mids <- dist - (incr/2)
  
  # Evaluate the lognormal function at the midpoint of each increment,
  # avoiding overflow for large values of the exponent
  es <- 0.5 * (log(mids/X0)/Xb)^2
  temp <- ifelse(es > 50, 0, exp(-es))
  
  # Convert to m2 by calculating area of the ring for each increment
  # To calculate the area of the ring, subtract the area of the circle
  # of increment i-1 from the circle of increment i
  temp = temp * pi * (dist^2 - (c(0, dist[1:(length(dist)-1)]))^2)
  
  # Sum it to get normalizer
  norm <- sum(temp)
  
  #########################
  # Calculate model predicted
  #########################
  predicted <- rowSums(TLP/norm * (dbhs/30.0)^alpha *
                         exp(-0.5*(log(distances/X0)/Xb)^2), na.rm=T)
  
  
  # Scale expected to size of the basket
  predicted <- predicted * basket.size
  
  # Expected can = 0, even with canopy trees in neighborhood, under the
  # lognormal function when X0 is large and Xb is very small. This gives
  # an underflow that should be set to a very small number so that the
  # likelihood calculation can be done...
  ifelse(predicted == 0, 0.0000001, predicted)
  
}


