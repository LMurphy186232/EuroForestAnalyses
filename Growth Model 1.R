###############################################################################
# Growth Model 1
# 
# For growth between 2002 and 2006:
# g = gmax * size * competition effect
# 
# where g is the actual growth, and gmax is the potential growth, NCI captures
# the relative social status (subordinate, co-dominant, dominant) and species
# composition of the relevant neighborhood. 
# 
# For the size effect, we'll use the log-normal model (for now), so
# size = gm * exp(-1/2) * [(ln(dbh/Xo)/Xb)]^2 
#
# where Xo is the dbh at which max. growth occurs, and Xb determines the 
# breadth (variability) of the function.
#
# For competition, we use their Neighborhood Competition Index (NCI)
# 
# competition effect = exp(-C* NCI^D)
#
# NCI(focal, k) = dbh(target, k) * sum(across species) * 
# sum(across number of neighbors within R) * lambda(species(focal), 
# species(competitor)) * eta(species k) * [dbh(neighbor, species)^alpha(species)
# /distance(neighbor, species)^beta(species)
# Eta(k) is a vector with the three social status designations. Similarly, we
# can assume that alpha and beta are the same in each species.
#
# Edit history:
# 1/6/2020 Lora Murphy - first draft
###############################################################################

library(likelihood)
projdir <- "C:/Users/lora/Documents/SORTIE/PROJECTS/AustriaPoland"

#-----------------------------------------------------------------------------#
# Input dataset. Fields:
# uid - unique identifier
# plot - study plot
# spid - species id (101: P. abies; 201: Abies alba; 1001: Fagus sylvatica)
# ht.2002: height estimate (in meters) 2002
# ht.2006: height estimate 2006
# soc2002: social position 2002 (d - dominant, cd - codominant, s - subordinate, 
# NA - unclear)
# soc2006: social position 2006
# x - in UTM (however, x-values are negative, not sure why. They map out roughly
# in the plot shape)
# y - in UTM
#-----------------------------------------------------------------------------#
setwd(projdir)
#dat <- read.csv("rothwald_joplot2.csv", as.is = T)
dat <- read.csv("rothwald_trees.csv")

# It is probable that we will want DBH in units of meters for NCI, so:
dat$dbh.m <- dat$dbh.2002 / 100.0

# Add diameter growth in mm/yr - this is making the assumption that all trees
# were censused at the same time for both censuses, and they were 4 years
# apart exactly. If this improbable situation is in fact not true, this is
# easily adjusted
dat$growth_mm_yr <- ((dat$dbh.2006 - dat$dbh.2002) * 10.0) / 4.0

# I'm going to drop a large outlier with 153 mm of radial growth
dat <- subset(dat, growth_mm_yr < 150)



#-----------------------------------------------------------------------------#
# Possible desired controls - currently set to not really remove trees. Just
# in case it makes sense to do further filtering here
#-----------------------------------------------------------------------------#
# Minimum target tree DBH, cm
min_target_dbh <- 0

# Minimum neighbor tree DBH, cm
min_neighbor_dbh <- 0

# Max radius to search for neighbors, m. 
##############################################
# Charlie suggested that this should become a variable that the model fits.
# A good suggestion for a future model iteration. We probably need to look at
# a few different values here and check the fitted results for alpha and beta
# in NCI to see if it looks like we've captured all the possible competitive
# neighbors.
##############################################
max_neighbor_radius <- 20

# Min growth threshold to be included, mm/yr
min_growth_allowed <- -25



#-----------------------------------------------------------------------------#
# It's easiest to analyze species one at a time, rather than trying to keep
# track of multiple values of each parameter. We really only have enough 
# individuals of species 1001 (Fagus sylvatica) to analyze here, so set that 
# as our analysis species. If we had multiple analyzable species in the dataset
# we could do something here to manage that.
#-----------------------------------------------------------------------------#
# Species we're doing right now
spp <- 1001




#-----------------------------------------------------------------------------#
# Subset to our target trees that we will analyze for growth. Targets are of
# the right species, right size, and don't have out-of-bounds growth values.
# This dataset is pretty clean so this isn't very meaningful here, but is an
# example of what we might need to do. Other conditions can also be added as
# needed.
#-----------------------------------------------------------------------------#
x <- which(dat$spid         == spp               & 
           dat$dbh.2002     >  min_target_dbh    &
           dat$growth_mm_yr >= min_growth_allowed)
targets <- dat[x,]





#-----------------------------------------------------------------------------#
# Get our set of potential neighbors, with whatever filtering is required 
# there.
#-----------------------------------------------------------------------------#
x <- which(dat$dbh.2002 > min_neighbor_dbh &
           !is.na(dat$soc2002))
neighbors <- dat[x,]


#-----------------------------------------------------------------------------#
# Now we need to collect information on the set of neighbors for each of our
# target trees. In our case, we need to know each neighbor's species, size, 
# distance, and social dstatus. 
# 
# To hold the neighbor info, prepare matrixes for each type of neighbor 
# information. Each row will be the neighbors of a single target, with the rows
# the same as in "targets". Each column is a neighbor. NAs are fully 
# acceptable as we expect different numbers of neighbors for different targets.
# 
# DBH is going to be in meters to make NCI units pleasant.
#-----------------------------------------------------------------------------#
# Function to calculate distance
neighdist<-function(targetx, targety, neighborx, neighbory) {
    sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)
}


# Declare the matrixes. We're sizing them as though every single tree is a 
# potential neighbor, which would be stupid on a larger dataset. Here, though,
# it lets us proceed without knowing how many neighbors we need to make room
# for, and then we are free to resize them later.
distances <- matrix(NA, nrow=nrow(targets), ncol=nrow(neighbors))
dbhs      <- matrix(NA, nrow=nrow(targets), ncol=nrow(neighbors))
species   <- matrix(NA, nrow=nrow(targets), ncol=nrow(neighbors))
status    <- matrix(NA, nrow=nrow(targets), ncol=nrow(neighbors))

# Save status as a factor. Now we can use the status value as an array index.
neighbors$soc2002 <- as.factor(neighbors$soc2002)
soc_levels <- levels(neighbors$soc2002)

# Save species as a factor too, same reason
neighbors$spid <- as.factor(neighbors$spid)
spp_levels <- levels(neighbors$spid)

# Populate the matrixes
for (i in 1:nrow(targets)) {
  
  # Only take trees within the neighborhood distance
  one_distances  <- neighdist(targets$x[i], targets$y[i], neighbors$x, neighbors$y)
  x <- which(one_distances <= max_neighbor_radius & one_distances > 0)
  
  distances[i,1:length(x)] <- one_distances    [x]
  dbhs     [i,1:length(x)] <- neighbors$dbh.m  [x]
  species  [i,1:length(x)] <- neighbors$spid   [x]
  status   [i,1:length(x)] <- neighbors$soc2002[x]
}

# Now we can resize the matrixes
x <- colSums(!is.na(distances))
distances <- distances[,(x>0)]
dbhs      <- dbhs     [,(x>0)]
species   <- species  [,(x>0)]
status    <- status   [,(x>0)]




#-----------------------------------------------------------------------------#
# Define the model
# Units of maxgrowth are mm/yr
#-----------------------------------------------------------------------------#
model1 <- function(maxgrowth, sizeX0, sizeXb, lambda, alpha, beta, 
                   eta, C, D)
{
    
  # Size effect. In this case, we can leave DBH in units of cm if we want
  size.effect <- exp(-0.5*(log(targets$dbh.2002/sizeX0)/sizeXb)^2)
    
  # Use the species as indexes to assemble a matrix of the correct lambda 
  # values
  lambda_vals <- lambda[species]
  dim(lambda_vals) <- dim(species)
  
  # Use the social status to assemble a matrix of the correct eta values
  eta_vals <- eta[status]
  dim(eta_vals) <- dim(status)
  
  # Now calculate NCI. Both distances and DBH are in meters here
  nci <- rowSums((lambda_vals * eta_vals * (dbhs ^ alpha)/(distances ^ beta)), na.rm=T)
  
  # Competition effect. Values of C are often very close to zero so we do this
  # to improve our ability to manage search ranges.
  competition_effect <- exp(-(C/100) * nci^D)
  
  # Final growth rate
  maxgrowth * size.effect * competition_effect
}



#-----------------------------------------------------------------------------#
# Analysis
#-----------------------------------------------------------------------------#

# Initial parameter values. In theory, these shouldn't matter too much as long
# as they are physically possible and within the search ranges. In practice, if
# it is difficult to make the model converge, one thing to try is playing 
# around with these values to make them actually reasonable. This helps in the
# case of a search surface that is covered with lots of pits that the search 
# gets trapped in. You know this happens when it keeps finding different, 
# terrible answers, all found pretty quickly
par <- list(
  maxgrowth = 10,
  sizeX0    = 30,
  sizeXb    = 2,
  C         = 100,
  D         = 1,
  alpha     = 2,
  beta      = 2,
  lambda    = rep(1, length(spp_levels)),
  eta       = rep(1, length(soc_levels)),
  sd        = 2) # sd is the standard deviation of the normal PDF, which we 
                 # estimate as a parameter rather than trying to specify it

# Upper and lower bounds. If any of the results bump up against these limits,
# expand them
par_hi <- list(
  maxgrowth = 500,
  sizeX0    = 500,
  sizeXb    = 50,
  C         = 100000,
  D         = 5,
  alpha     = 4,
  beta      = 4,
  lambda    = rep(1, length(spp_levels)),
  eta       = rep(1, length(soc_levels)),
  sd        = 1000)

par_lo <- list(
  maxgrowth = 0,
  sizeX0    = 0,
  sizeXb    = 0.2,
  C         = 0,
  D         = 1,
  alpha     = 0,
  beta      = 0,
  lambda    = rep(0, length(spp_levels)),
  eta       = rep(0, length(soc_levels)),
  sd        = 0.2)


# Assign some names to make the values for lambda and eta easy to identify
names(par$lambda) <- spp_levels
names(par$eta) <- soc_levels

# This is harder to explain but this is the list of parameter values that we 
# are not fitting, so they are separated out. These are all used by dnorm
var <- list(mean = "predicted", x = "growth_mm_yr", log=T)

# The annealing schedule (initial temperature and rate of reduction) is
# not fixed and can be played with. Hotter initial temperatures cause a
# more wide-ranging initial search that takes longer to settle down, 
# which can be good if it's a rather complicated function. A slower rate
# of temperature drop keeps the search going longer. You just want to
# be sure, looking at the annealing history, that the search had enough
# "energy" to have found a good answer and long enough to have zeroed in
# on it. I usually want at least 2 temperature drops without a difference
# in likelihood to be sure that the answer was found.
results<-anneal(model = model1, par = par, var = var, source_data = targets, 
                par_lo = par_lo, par_hi = par_hi, pdf = dnorm, 
                dep_var = "growth_mm_yr", initial_temp = 6, temp_red = 0.95,
                max_iter = 5000)

par <- results$best_pars
res <- support_limits(model = model1, par = par, var = var, source_data = targets, 
                       par_lo = par_lo, par_hi = par_hi, pdf = dnorm)
all_good <- length(res) == 2
while (!all_good) {
  par[[res$i]][res$j] <- res$parval
  res <- support_limits(model = model1, par = par, var = var, source_data = targets, 
                        par_lo = par_lo, par_hi = par_hi, pdf = dnorm)
  all_good <- length(res) == 2
}  
results<-anneal(model = model1, par = par, var = var, source_data = targets, 
                par_lo = par_lo, par_hi = par_hi, pdf = dnorm, 
                dep_var = "growth_mm_yr", initial_temp = 6, temp_red = 0.95,
                max_iter = 500)
res <- support_limits(model = model1, par = par, var = var, source_data = targets, 
                      par_lo = par_lo, par_hi = par_hi, pdf = dnorm)
results$best_pars <- par
results$upper_limits <- res$upper_limits
results$lower_limits <- res$lower_limits


write_results(results, "Growth Model 1 R results.txt",data=F, print_whole_hist=F)

save(results,file="Growth Model 1 R results.Rdata")




#-----------------------------------------------------------------------------#
# Plotting
#-----------------------------------------------------------------------------#
do_post <- F
if (do_post) {
  windows()
  plot(results$source_data$growth_mm_yr, results$source_data$predicted)
  
  # Look at the size effect predicted
  x <- 1:100
  size.effect <- exp(-0.5*(log(x/results$best_pars$sizeX0)/results$best_pars$sizeXb)^2)
  plot(x, size.effect,
       type="l", lwd=3,
       main = "Size effect for Mayer-Wegelin, model 1",
       xlab = "2002 DBH, cm")
  
  # Check dbh vs growth, laid on top
  # Relativize growth so it fits on top
  points(dat$dbh.2002, dat$growth_mm_yr/max(dat$growth_mm_yr), pch=20)
  
  # What kind of range are we getting on NCI values?
  lambda_vals <- results$best_pars$lambda[species]
  dim(lambda_vals) <- dim(species)
  
  # Use the social status to assemble a matrix of the correct eta values
  eta_vals <- results$best_pars$eta[status]
  dim(eta_vals) <- dim(status)
  
  # Now calculate NCI. Both distances and DBH are in meters here
  nci <- rowSums((lambda_vals * eta_vals * (dbhs ^ results$best_pars$alpha)/(distances ^ results$best_pars$beta)), na.rm=T)
  
  # Competition effect. Values of C are often very close to zero so we do this
  # to improve our ability to manage search ranges.
  competition_effect <- exp(-(C/100) * nci^D) 
  
}
