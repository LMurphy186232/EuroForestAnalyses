###############################################################################
# Disperse shell
# 
# Filling out the pieces needed for the disperse model...
# 
# Lora Murphy, 5/20/2020
###############################################################################

library(likelihood)
projdir <- "C:/Users/lora/Documents/SORTIE/PROJECTS/AustriaPoland"

#-----------------------------------------------------------------------------#
# Controls. Could also become fitted variables
#-----------------------------------------------------------------------------#
# Minimum parent tree DBH, cm
min_parent_dbh <- 0

# Max radius to search for parents, m. 
max_parent_radius <- 20


#----- Controls for the normalizer -------------------------------------------#
# Width of rings, in whatever distance units we're using
norm_incr <- 0.25
# Length of the entire disperse kernel, in whatever distance units. This is
# more about capturing the long tail of the function, to approximate a true
# integral, than the distance that we actually think trees will disperse
norm_maxdist <- 250
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# Data prep
#-----------------------------------------------------------------------------#





#-----------------------------------------------------------------------------#

#----- Controls for the normalizer -------------------------------------------#
# Width of rings, in whatever distance units we're using
incr <- 0.25
# Length of the entire disperse kernel, in whatever distance units
maxdist <- 250


#-----------------------------------------------------------------------------#
# Function for an isotropic lognormal normalizer.
# Lognormal: exp(-0.5(ln(dist/X0)/Xb))
# The procedure is to divide a circle into a series of rings of constant width
# (size of your choice). At the distance of the midpoint of each ring, the 
# lognormal function is calculated, and the result multiplied by the area of 
# the ring. The final result is the sum of these values across all rings.
#-----------------------------------------------------------------------------#
lognormal_iso_normalizer <- function(X0, Xb) {

  # Get the distances to each ring
  dist <- seq(from=norm_incr, to=norm_maxdist, by=norm_incr)
  
  # Get the midpoint of each ring - halfway between the increments
  mids <- dist - (norm_incr/2)
  
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

#-----------------------------------------------------------------------------#
# Function for an isotropic weibull normalizer.
###############################
# NOTE: Charlie had some extra terms in his Weibull that I don't understand.
# They may be some sort of scaling thing. We'll have to get his input.
###############################

# Weibull: exp((-D*dist^theta))
# Arguments:
# D, theta: parameters for weibull.
#-----------------------------------------------------------------------------#
weibull_iso_normalizer <- function(D, theta) {
  
  # Get the distances to each ring
  dist <- seq(from=norm_incr, to=norm_maxdist, by=norm_incr)
  
  # Get the midpoint of each ring - halfway between the increments
  mids <- dist - (norm_incr/2)

  # Evaluate the weibull function at the midpoint of each increment
  es <- D * dist^theta
  
  # Avoid overflow for large values of the exponent
  temp <- ifelse(es > 50, 0, exp(-es))
  
  # Convert to units^2 by calculating area of the ring for each increment
  # To calculate the area of the ring, subtract the area of the circle
  # of increment i-1 from the circle of increment i
  temp = temp * pi * (dist^2 - (c(0, dist[1:(length(dist)-1)]))^2)
  
  # Sum it to get normalizer
  sum(temp)
  
}


#-----------------------------------------------------------------------------#
# Lognormal isotropic model
#-----------------------------------------------------------------------------#
lognormal_iso_model <- function(TSP, alpha, X0, Xb) {
  
  # Calculate normalizer
  norm <- lognormal_iso_normalizer(X0, Xb)
  
  predicted <- rowSums(TSP/norm * (dbhs/30.0)^alpha *
                         exp(-0.5*(log(distances/X0)/Xb)^2), na.rm=T)
  
  
  # Expected can = 0, even with canopy trees in neighborhood, under the
  # lognormal function when X0 is large and Xb is very small. This gives
  # an underflow that should be set to a very small number so that the
  # likelihood calculation can be done...
  ifelse(predicted == 0, 0.0000001, predicted)
  
}
