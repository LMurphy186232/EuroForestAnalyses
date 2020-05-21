###############################################################################
# Disperse shell
# 
# Filling out the pieces needed for the disperse model...
# 
# Lora Murphy, 5/20/2020
###############################################################################

incr <- 0.25
maxdist <- 250

#-----------------------------------------------------------------------------#
# Function for an isotropic lognormal normalizer.
# The procedure is to divide a circle into a series of rings of constant width
# (size of your choice). At the distance of the midpoint of each ring, the 
# lognormal function is calculated, and the result multiplied by the area of 
# the ring. The final result is the sum of these values across all rings.
# 
# Arguments:
# incr: width of the rings.
# maxdist: radius of circle = maximum possible disperse distance.
# X0, Xb: parameters for lognormal.
#-----------------------------------------------------------------------------#
lognormal_iso_normalizer <- function(incr, maxdist, TLP, alpha, X0, Xb) {
  
  
  # Get the distances to each ring
  dist <- seq(from=incr, to=maxdist, by=incr)
  
  # Get the midpoint of each ring - halfway between the increments
  mids <- dist - (incr/2)
  
  # Evaluate the lognormal function at the midpoint of each increment
  es <- 0.5 * (log(mids/X0)/Xb)^2
  
  # Avoid overflow for large values of the exponent
  temp <- ifelse(es > 50, 0, exp(-es))
  
  # Convert to units^2 by calculating area of the ring for each increment
  # To calculate the area of the ring, subtract the area of the circle
  # of increment i-1 from the circle of increment i
  temp = temp * pi * (dist^2 - (c(0, dist[1:(length(dist)-1)]))^2)
  
  # Sum it to get normalizer
  sum(temp)
  
}

