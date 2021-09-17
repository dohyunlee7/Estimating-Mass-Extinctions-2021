

# read in the functions that calculate the CI
# change pathname to match your file location
source("~/Desktop/duration28c.R")


# read in dataset of Seymour Island ammonites
# change pathname to match your file location
data <- matrix(scan("~/Desktop/seymourNA.txt"), 24, 10, byrow=T)  
data[data<1000] <- NA       # use only values above 1000 m
data <- data-1000           # rescale data so base = 0


# run the algorithm
par(mfrow=c(1,1), omi=c(1,1,1,1)*.2)
deltaCI28c(data, stepsize=4, PLOT=1, xplotmax=67, yplotmax=270)


# the function should return "0.0 62.0 57.5", or something close to that
#  (subject to random simulation variability)
# this means the CI for delta is (0, 62), and the observed value of d is 57.5

