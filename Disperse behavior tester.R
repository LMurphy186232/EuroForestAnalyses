###############################################################################
# Test the approach that we're outlining for the disperse behavior. 
# 
# For a tree map, this will create a time series of seed numbers so they can
# be checked for reasonableness. Obviously all other tree lifecycle dynamics
# are ignored.
###############################################################################

#----- Plot data -------------------------------------------------------------#
# Read in some tree data - really only to get a realistic DBH distribution.
dat <- read.csv("rothwald_joplot.csv")
dat <- dat[dat$spid == 1001, ]

# Planned mean STR time series - one value per year
mean_STR <- runif(n = 100, min = 1000, max = 10000)

# Proportion reproducing time series - one value per year
# This is where we can add detail later
prop_reproducing <- runif(n = 100)
#-----------------------------------------------------------------------------#


#----- Tree level parameters -------------------------------------------------#
# Assign a mean rho value for each tree that will be used each year. How
# does this get assigned?
rho_mean <- runif(n = nrow(dat))

# Keep a memory of the STR because the previous year's value is used in the
# current year. 
str <- matrix(0, nrow = nrow(dat), ncol = length(mean_STR) + 1)
colnames(str) <- paste0("STR_t_", 0:length(mean_STR))

# STR initial conditions - year 0. HOW DO WE GET THE VALUE OF THE PREVIOUS 
# YEAR FOR THE FIRST YEAR OF THE TIME SERIES? Just using pop mean right now
str[,1] <- mean_STR[1]
#-----------------------------------------------------------------------------#



#----- Matrixes to store values ----------------------------------------------#
# Keep a record of rho for each year to examine later and aid in 
# troubleshooting - SORTIE won't keep this
rho <- matrix(0, nrow = nrow(dat), ncol = length(mean_STR))
colnames(rho) <- paste0("rho_t_", 1:length(mean_STR))

# Here's where we'll put the seed count
seed_counts <- matrix(0, nrow = nrow(dat), ncol = length(mean_STR))
colnames(seed_counts) <- paste0("seeds_t_", 1:length(mean_STR))
#-----------------------------------------------------------------------------#



#----- Population level parameters -------------------------------------------#
# Omega, the population-level scaling factor. This needs either a value or a 
# rule for calculating.
omega <- 1

# Beta, the size parameter. 
beta <- 1

#-----------------------------------------------------------------------------#
# Do the simulation
#-----------------------------------------------------------------------------#
for (yr in 1:length(mean_STR)) {
  
  #---------------------------------------------------------------------------#
  # Get each tree's STR value for this year
  #---------------------------------------------------------------------------#
  
  # Draw the rho value for this year. What is the distribution? What are the
  # other parameters of the distribution (in this case, standard deviation)
  # and where do they come from?
  # 
  # Enforcing to be strictly positive with a definitely uncool operation here.
  # Will have to figure out what it means if rho is negative, if drawn from a
  # distributiom where it can be negative
  rho[,yr] <- abs(rnorm(n = nrow(rho), mean = rho_mean, sd = 1))
  
  # Calculate STR. Remember that for the STR matrix, we are also preserving
  # year zero, so there is one extra column
  str[,(yr+1)] <- mean_STR[yr] + (rho[,yr] * str[,yr])
  #---------------------------------------------------------------------------#
  
  
  
  #---------------------------------------------------------------------------#
  # Calculate final number of seeds
  #---------------------------------------------------------------------------#
  seed_counts[,yr] <- str[,(yr+1)] * (dat$dbh.2006^beta) * omega
  
  # Set those not participating to zero. This isn't how SORTIE will go about
  # this but it's functionally equivalent
  seed_counts[,yr] <- ifelse(runif(n = nrow(seed_counts)) < prop_reproducing[yr],
                             0, seed_counts[yr,])
  #---------------------------------------------------------------------------#
}
#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
# Quick graphs here
#-----------------------------------------------------------------------------#
# Actual proportion reproducing
actual_proportion <- apply(seed_counts, MARGIN=2, FUN = function(x){sum(x > 0)/length(x)})
plot(1:length(actual_proportion), actual_proportion, type="l",
     lwd = 3, col = "black",
     ylim = c(0,1),
     main = "Proportion of trees reproducing",
     xlab = "Year", ylab = "Proportion reproducing")
lines(1:length(prop_reproducing), prop_reproducing,
      lwd = 1, col="darkgray")
legend("topright", legend = c("Actual", "Requested"),
       lwd = c(3, 1), col=c("black", "darkgray"))