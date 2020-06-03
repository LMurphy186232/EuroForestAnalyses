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

