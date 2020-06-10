###############################################################################
# Functions for disperse
# 
# Lora Murphy, 6/3/2020
###############################################################################


#-----------------------------------------------------------------------------#
# Zero-inflated normal PDF
# zprob is the probability of a zero
# varm is a parameter that controls the standard deviation as a linear functio
# Normally we would use a discrete PDF for seed data like the poisson or neg
# binomial, but we are using the seeds per m2 here since the normalizer is in
# units of square meters. A discrete could be used along with seed numbers if 
# the normalizer is adjusted based on trap size.
#-----------------------------------------------------------------------------#
zinf_dnorm <- function(zprob, mean, varm)
{ 
  sd <- mean * varm
  loglh <- log(ifelse(seeds$seeds_m2 == 0,
              zprob + (1-zprob) * dnorm(0,              mean, sd, log=F), 
                      (1-zprob) * dnorm(seeds$seeds_m2, mean, sd, log=F)))
  return(loglh) 
}
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# Function for an isotropic lognormal normalizer.
# Lognormal: exp(-0.5(ln(dist/X0)/Xb))
#
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





#-----------------------------------------------------------------------------#
# Zero-inflated heteroscedastic normal PDF
# zprob is the probability of a zero
# sigma is a parameter that controls the standard deviation as a linear function
# of the mean
# Normally we would use a discrete PDF for seed data like the poisson or neg
# binomial, but we are using the seeds per m2 here since the normalizer is in
# units of square meters. A discrete could be used along with seed numbers if 
# the normalizer is adjusted based on trap size.
#-----------------------------------------------------------------------------#
zinf_dnorm_v <- function(zprob, mean, sigma)
{ 
  sd <- mean * sigma
  loglh <- log(ifelse(seeds$seeds_m2 == 0,
              zprob + (1-zprob) * dnorm(0,              mean, sd, log=F), 
                      (1-zprob) * dnorm(seeds$seeds_m2, mean, sd, log=F)))
  return(loglh) 
}
#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
# Zero-inflated heteroscedastic normal PDF, zprob is by year
# zprob is the probability of a zero
# sigma is a parameter that controls the standard deviation as a linear function
# of the mean
# Normally we would use a discrete PDF for seed data like the poisson or neg
# binomial, but we are using the seeds per m2 here since the normalizer is in
# units of square meters. A discrete could be used along with seed numbers if 
# the normalizer is adjusted based on trap size.
#-----------------------------------------------------------------------------#
zinf_dnorm_v_yr <- function(zprob, mean, sigma)
{ 
  sd <- mean * sigma
  loglh <- log(ifelse(seeds$seeds_m2 == 0,
      zprob[seeds$yr] + (1-zprob[seeds$yr]) * dnorm(0,              mean, sd, log=F), 
                        (1-zprob[seeds$yr]) * dnorm(seeds$seeds_m2, mean, sd, log=F)))
  return(loglh) 
}
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# Function to calculate the area of intersection of circle and square. Got 
# this off of Ars Technica because I didn't feel like working out the formula 
# on my own. The basic approach here is to consider each quadrant of the circle
# separately and sum the four areas.
# 
# I didn't work out all the math independently but I wrote a simple numerical
# approximator and compared a bunch of test cases (see separate testing script)
# and this gives the correct answers. The approximator would work too but this
# gives the exact right answer, which is more awesomer.
# 
# P_rect_top_left(x1, y1), P_rect_bottom_right(x2, y2), 
# P_circle_center(mx, my) and radius(r).
#-----------------------------------------------------------------------------#
frac_in_square <- function(x1, y1, x2, y2, mx, my, r) {
  
  x1 <- x1 - mx
  x2 <- x2 - mx
  y1 <- y1 - my
  y2 <- y2 - my
  
  # This does each quadrant separately. Original script didn't use the abs() but I found
  # it was necessary.
  a <- abs(sA(r, x2, y1) - sA(r, x1, y1) - sA(r, x2, y2) + sA(r, x1, y2))
  
  # Calculate the ratio of the area of intersection to the area of the whole circle
  a/(pi*r^2)
}

sA <- function(r, x, y) {
  
  if (x < 0) {
    return(-sA(r, -x, y))
  }
  
  if (y < 0) {
    return(-sA(r, x, -y))
  }
  
  if (x > r) {
    x = r
  }
  
  if (y > r) {
    y = r
  }
  
  # If square is bigger than circle
  if (x*x + y*y > r*r) {
    a = r*r*asin(x/r) + x*sqrt(r*r-x*x) + 
      r*r*asin(y/r) + y*sqrt(r*r-y*y) - 
      r*r*pi/2
    
    a <- a * 0.5
    
  } else {
    # Circle is bigger than square
    a = x*y
  }
  a
}

