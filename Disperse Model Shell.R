###############################################################################
# Disperse shell
# 
# Filling out the pieces needed for the disperse model.
# 
# Issues to address: 
# Size of potential parents. We have DBH measurements in 2002, 2006 (with some
# missing), and 2012. We have seed data from 2003 - 2018. Do you grab the 
# nearest year? The most recent census? Do you interpolate missing DBH values
# for trees that you know to be alive in 2006? I am going to do a default of
# grabbing the most recent previous census with a value (going back to 2002
# if 2006 is missing). This is just to get something started. It isn't a 
# recommendation on my part.
# 
# I am throwing all the years of seed data together.
# 
#
# We could separate the functions into a different file and source them if
# desired.
# 
# Lora Murphy, 5/20/2020
###############################################################################

library(likelihood)

# Obviously the directories reflect how I set this up on my own machine so
# adjust as necessary
projdir <- "C:/Users/lora/Documents/SORTIE/PROJECTS/AustriaPoland/"

# Pull in our functions
source("disperse_functions.R")

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
# Data prep - seeds
# Species codes:
# Bu - beech Fagus sylvatica
# Ta - silver fir Abies Albany
# Fi - Spruce Picea abies
# Os - other species
# 
# We appear to be missing location info for seed trap MW27
#-----------------------------------------------------------------------------#
# Seed data
setwd(paste(projdir, "Disperse", sep="/"))
seeds <- read.csv("seeds_by_trap_year.csv", as.is=T)
# Get rid of what appear to be row numbers
seeds$X <- NULL

# Mayer-Wegelin only
x <- grep("mw", seeds$uid)
seeds <- seeds[x,]

# Which species do I want to do? In this case beech
seeds <- seeds[seeds$species == "BU",]

# Add locations
trap_loc <- read.csv("MW_seedtrap_coords.csv")
trap_loc$uid <- paste0("mw", trap_loc$id_trap)
ind <- match(seeds$uid, trap_loc$uid)
seeds$x <- trap_loc$x[ind]
seeds$y <- trap_loc$y[ind]

# No NAs
x <- which(!is.na(seeds$seeds_abs) & !is.na(seeds$x))
seeds <- seeds[x,]
#-----------------------------------------------------------------------------#





#-----------------------------------------------------------------------------#
# Data prep - adults
# spid - species id (101: P. abies; 201: Abies alba; 1001: Fagus sylvatica)
# 
# I am uncomfortable with the number of NAs in the alive.2012 field, and the
# missing 2006 records. We will definitely be missing some parents here.
#-----------------------------------------------------------------------------#
setwd(projdir)
adults <- read.csv("rothwald_trees.csv")

# There are trees marked as dead in an early census and alive in a later one.
# Correct to assume no resurrections actually occurred
x <- which(adults$alive.2012)
adults$alive.2002[x] <- T
adults$alive.2006[x] <- T

x <- which(adults$alive.2006)
adults$alive.2002[x] <- T


# We are missing a lot of records from 2006. Each seed observation is going to
# get the most recent previous DBH measurement for parents. So fill in missing
# DBH values from a previous census if available. 
adults$dbh.2006 <- ifelse(is.na(adults$dbh.2006), adults$dbh.2002, adults$dbh.2006)
adults$dbh.2012 <- ifelse(is.na(adults$dbh.2012), adults$dbh.2006, adults$dbh.2012)

# Only Mayer-Wegelin beeches
adults <- adults[adults$spid == 1001 & adults$plot == "mw",]

# Construct 3 parent datasets: for 2002, 2006, and 2012. Many different 
# possible approaches here for making it easy to figure out which trees are
# parents.
adults_by_year = list()

x <- which(!is.na(adults$dbh.2002 & adults$alive.2002 == T))
adults_by_year[["2002"]] <- data.frame(DBH = adults$dbh.2002[x], 
                                         X = adults$x_loc   [x], 
                                         Y = adults$y_loc   [x])

x <- which(!is.na(adults$dbh.2006 & adults$alive.2006 == T))
adults_by_year[["2006"]] <- data.frame(DBH = adults$dbh.2006[x], 
                                         X = adults$x_loc   [x], 
                                         Y = adults$y_loc   [x])

x <- which(!is.na(adults$dbh.2012 & adults$alive.2012 == T))
adults_by_year[["2012"]] <- data.frame(DBH = adults$dbh.2012[x], 
                                         X = adults$x_loc   [x], 
                                         Y = adults$y_loc   [x])
#-----------------------------------------------------------------------------#









#-----------------------------------------------------------------------------#
# Construct the master dataset. Each trap, for each year, is its own 
# observation. For each, identify potential parents and their distance and most
# recent size. This is the same approach as for growth models with NCI 
# competition.
#-----------------------------------------------------------------------------#
# Function to calculate distance
neighdist<-function(targetx, targety, neighborx, neighbory) {
  sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)
}



# Declare the matrixes. We're sizing them as though every single tree is a 
# potential parent, which would be stupid on a larger dataset. Here, though,
# it lets us proceed without knowing how many neighbors we need to make room
# for, and then we are free to resize them later.
distances <- matrix(NA, nrow=nrow(seeds), ncol=nrow(adults))
dbhs      <- matrix(NA, nrow=nrow(seeds), ncol=nrow(adults))

# Populate the matrixes
for (i in 1:nrow(seeds)) {
  
  year = "2012"
  if (seeds$year[i] <  2006)                         year = "2002"
  if (seeds$year[i] >= 2006 && seeds$year[i] < 2012) year = "2006"
  
  # Only take trees within the neighborhood distance
  one_distances  <- neighdist(seeds$x[i], seeds$y[i], 
                              adults_by_year[[year]]$X, 
                              adults_by_year[[year]]$Y)
  x <- which(one_distances <= max_parent_radius & one_distances > 0)
  
  distances[i,1:length(x)] <- one_distances             [x]
  dbhs     [i,1:length(x)] <- adults_by_year[[year]]$DBH[x]
}

# Now we can resize the matrixes
x <- colSums(!is.na(distances))
distances <- distances[,(x>0)]
dbhs      <- dbhs     [,(x>0)]
#-----------------------------------------------------------------------------#





#-----------------------------------------------------------------------------#
# Possible correction approach for traps near the plot edge that might not get
# to search all possible parents: there is a function at the end of this page
# that will calculate the fraction of your neighborhood found within the plot.
# Tuck this into the above code and do something with frac.
#-----------------------------------------------------------------------------#
frac <- frac_in_square(0, 0, 100, 100, seeds$x[i], seeds$y[i], max_parent_radius)
#-----------------------------------------------------------------------------#





#----- Controls for the normalizer -------------------------------------------#
# Width of rings, in whatever distance units we're using
incr <- 0.25
# Length of the entire disperse kernel, in whatever distance units
maxdist <- 250














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
#-----------------------------------------------------------------------------#






#-----------------------------------------------------------------------------#
# Weibull isotropic model
#-----------------------------------------------------------------------------#
weibull_iso_model <- function(TSP, alpha, D, theta) {
  
  # Calculate normalizer
  norm <- weibull_iso_normalizer(D, theta)
  
  predicted <- rowSums(TSP/norm * (dbhs/30.0)^alpha *
                         exp(-D*distances^theta), na.rm=T)
  
  
  # Guard against underflows
  ifelse(predicted == 0, 0.0000001, predicted)
  
}
#-----------------------------------------------------------------------------#






#-----------------------------------------------------------------------------#
# Weibull isotropic model
#-----------------------------------------------------------------------------#
weibull_iso_model <- function(TSP, alpha, D, theta) {
  
  # Calculate normalizer
  norm <- weibull_iso_normalizer(D, theta)
  
  predicted <- rowSums(TSP/norm * (dbhs/30.0)^alpha *
                         exp(-D*distances^theta), na.rm=T)
  
  
  # Guard against underflows
  ifelse(predicted == 0, 0.0000001, predicted)
  
}
#-----------------------------------------------------------------------------#






#-----------------------------------------------------------------------------#
# Weibull isotropic model, minimum size fitted as parameter; add minsize to
# par, par_lo, and par_hi
#-----------------------------------------------------------------------------#
weibull_iso_model_var_minsize <- function(TSP, alpha, D, theta, minsize) {
  
  # Calculate normalizer
  norm <- weibull_iso_normalizer(D, theta)
  
  # Set everything below min size to NA
  dbhtemp[[which(dbhs < minsize)] <- NA
  
  predicted <- rowSums(TSP/norm * (dbhtemp/30.0)^alpha *
                         exp(-D*distances^theta), na.rm=T)
  
  
  # Guard against underflows
  ifelse(predicted == 0, 0.0000001, predicted)
  
}
#-----------------------------------------------------------------------------#






#-----------------------------------------------------------------------------#
# Annealing
#-----------------------------------------------------------------------------#
par    <- list(TSP = 600,  alpha = 2,     X0 = 15,     Xb = 0.5, zprob = 0.5, sd = 2)
par_lo <- list(TSP = 0,    alpha = 0.001, X0 = 0.0001, Xb = 0.5, zprob = 0,   sd = 0.2)
par_hi <- list(TSP = 5000, alpha = 4,     X0 = 100,    Xb = 20,  zprob = 1,   sd = 1000)


var <- list(mean = "predicted")

windows()
results<-anneal(model = lognormal_iso_model, par = par, var = var, 
                source_data = seeds, par_lo = par_lo,
                par_hi = par_hi, pdf = zinf_dnorm, 
                dep_var = "seeds_m2", 
                initial_temp = 3, temp_red = 0.9, max_iter = 100000)
write_results(results, "results.txt")

#-----------------------------------------------------------------------------#






#-----------------------------------------------------------------------------#
# Setting up multiple values per variable: in this case, a zero probability
# for each year
#-----------------------------------------------------------------------------#
seeds$yr <- as.factor(seeds$year)
num_years <- length(levels(seeds$yr))

par    <- list(TSP = 600,  alpha = 2,     X0 = 15,     Xb = 0.5, zprob = rep(0.5, num_years), sd = 2)
par_lo <- list(TSP = 0,    alpha = 0.001, X0 = 0.0001, Xb = 0.5, zprob = rep(0,   num_years), sd = 0.2)
par_hi <- list(TSP = 5000, alpha = 4,     X0 = 100,    Xb = 20,  zprob = rep(1,   num_years),   sd = 1000)

# Add some names so 
names(par$zprob) <- paste0("zprob", as.character(levels(seeds$yr)))
var <- list(mean = "predicted")

windows()
results<-anneal(model = lognormal_iso_model, par = par, var = var, 
                source_data = seeds, par_lo = par_lo,
                par_hi = par_hi, pdf = zinf_dnorm_v_yr, 
                dep_var = "seeds_m2", 
                initial_temp = 3, temp_red = 0.9, max_iter = 100000)
write_results(results, "results.txt")

#-----------------------------------------------------------------------------#




