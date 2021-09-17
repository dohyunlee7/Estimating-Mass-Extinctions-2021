
# ---------------------------------------------------------------------------- #
#
# Confidence interval for extinction duration
# from Confidence intervals for the duration of a mass extinction
# by Steve Wang, Phil Everson, Brendan McVeigh, Aaron Zimmerman, Heidi Wong
# Paleobiology 38:2, 2012
# updated version: allows for non-uniform preservation
# additional code by Madison Shoraka, Melissa Zavez, Kevin Choi, Eric Zhang, 
#                    Danielle Rossetti Dos Santos, Peem Lerdputtipongporn
#
# ---------------------------------------------------------------------------- #



# --------------------------- Version history -------------------------------- #

# v4 6/10/08: added check for proposed extincton time < range extension
# v5 6/11/08: halts if all proposed extincton time > range extension
#           : also fixes problem if multiple taxa have same lowest y_i
# v6 6/13/08: converted into a function
# v7 6/14/08: added simulation loop
# v8 6/18/08: modified simulation for true thetas
# v9 10/20/08: added Phil’s joint posterior method
# v10 10/21/08: added constraint for max chi-square from any single taxon
# v11 2/27/09: took out second-order-stat constraint
# v13 6/12/10: added constraint for min chi-square
# v14 6/15/10: made brute force a function
# v15 7/6/10: new deltaCI simulating from delta scenarios
# v16 7/7/10: used theta-hat instead of y for d and delta; 1-tailed; from Phil
# v17 7/8/10: new way of choosing thetas for a given delta, from Brendan & Aaron
# v18 7/10/10: combines elements of v16 and v17 for setting hi/lo thetas; 
#              adds range extensions
# v19 7/19/10: added separate code for obs and obsdprime; 
#              fixed error in assigning intermediate thetas;
#              improved output        
# v20 7/26/10: went back to earlier method of assigning intermediate thetas        
# v21 7/27/10: added backtracking if obs d is not in the CI; added range charts
# v22 8/3/10: in backtracking, changed order of assigning thetahi and taxonhi
# v23 8/5/10: steps forward through all delta values instead of backtracking
# v24 8/6/10: samples thetas from beta distribution instead of uniform
# v25 8/7/10: determines thetahi/thetalo by adding to the highest/lowest finds
#             in order to avoid negative thetalo (instead of argmax thetahi/lo)
# v26 8/15/10: moved simulation code into a separate function ‘simdvalues’
# v27 8/22/10: in simtruehoriz, made MINTHETA/MAXTHETA the actual min and max
#              theta values, not just the min/max bounds for possible thetas
# v27a 12/16/10: final testing version of data-matched thetas
# v28 12/11/10: improved plot boundaries, added Meishan data (not data-matched)
# v29a 12/20/10: same as v28 but with old ‘efficient algorithm’ for testing (branch)
# v29b 12/29/10: checks percentage of total chi-sq attributable to one taxon
# v30 6/20/18: chooses lo theta and calculates hi theta as lo + delta (instead lo = hi - delta)
# v31 2/19/19: adds lower/upperlim on looping based on Madison and Melissa's idea
# v32 3/14/19: Madison and Melissa's binary search version (what they called 'No_Probabilities.R')
#              uses nsims1 and nsims2 for lower and upper endpoint
# v33 6/02/19: Eric and Kevin added non-uniform recovery
# v34 6/12/19: Eric and Kevin added lambda shrinking to non-uniform recovery
# v35 6/25/19: When choosing thetahi we use sum(n) instead of namx
# v35.1 6/25/19: When choosing thetahi we use the new PDF we derived
# v36 7/9/19: Estimate d proportional to ymax instead of estimating based on the value of d
# v37 7/9/19: Use the scaled d value for the lower while loop and the regular method for the upper bound
# v38 6/8/20: Cleaned up from v35.1. Standardized variable names, took out old code (e.g. accept),
#             keeps track of step #, improves plotting of binary search convergence, etc.
#             note: all v38 simulates d values for a given delta
# v38.1 6/08/20: starts at delta = 0; does not use min and max bounds
# v38.2 8/10/20: adds min and max search bounds back in as search limits (not as CI bounds)
# v38.3 3/03/21: uses min and max as search bounds, no simulation
# v38.4 3/04/21: uses min and for lower bound, simulation for upper bound
# v39 7/14/20: uses LRT instead of simulating d values given delta, by Peem and Danielle
# v40 goes back to v38: simulating d values instead of using LRT, by Emma and Jiung
# v40.1 4/09/21: uses signed version of d depending on correlation of ordering of taxa
# v40.2 4/23/21: multiplies d by signed correlation 
# v40.3 4/30/21: uses signed version of d depending on correlation but modifies rejection criteria
# v41.0 6/03/21: SCW implemented Fay et al 2015 method. Assumes only 2 taxa.
# v41.1 6/04/21: adjusted for >2 taxa by choosing max and min ci. Also adjusts CI to
#                work for abs(theta2-theta1) rather than (theta2 - theta1).
# v41.2 6/07/21: uses Fay for lower endpoint, simulation for upper
# v41.3 6/12/21: calculates Fay using normal and reversed order hi/lo, for debugging only
# v41.4 6/18/21: Fay method, revised and corrected




# ------------------------- Function definitions ------------------------------ #


# read in code for adaptive beta method CIs
#source("~/Documents/Research/Adaptive range extensions/ abm38.2.R")
source("~/Desktop/abm38.2.R")



# Draw range chart
plotrangechart <- function(data, ci, thetahi, thetalo)  {
  ord <- order(data[1,])
  data <- data[,ord]           # sort data by highest find
  ci <- ci[ord]
  ntaxa <- dim(data)[2]
  maxfinds <- dim(data)[1]
  ymax <- max(data[1,])
  plot( NULL, pch="", xlab="", ylab="", xlim=c(0, ntaxa+1), 
        ylim=c(-.3, max(ci)), bty="L", xaxt="n" )
  for(taxon in 1:ntaxa)  {
    # draw taxon lines
    lines( c(taxon,taxon), c(0,data[1,taxon]), lwd=.9, col="darkgray" )
    # plot fossil horizons
    points(rep(taxon,maxfinds), data[,taxon], pch=16, cex=.65)
    # plot confidence interval tops
    points(taxon, ci[taxon], pch="^", cex=.65, col="red")
    # plot theta hi values
    lines(taxon, thetahi, cex=.65, col = "green")
    #plot theta lo values
    #points(taxon, thetalo, cex=.65, col = "orange")
  }
}



# Simulate a sample from a reflected beta distribution with the given lambda
rrefbeta <- function(n, lambda)  {
  if(lambda<=0)  { 
    return(rbeta(n, 1, 1-lambda))
  }  else  return(rbeta(n, 1+lambda, 1))
}



# Simulate the gap between the highest fossil find and the true extinction time
simulategap <- function(lambda, n) {
  g <- function(x, lambda, n) (1-(1-x)^(1/n))^(1/(1-lambda))
  x <- runif(n=1)
  return(g(x, lambda, n))
}



# Simulate a dataset from thetas having true delta equal to ‘delta’, 
#   and calculate the d value from the dataset.
# Simulates only the highest finds, not the entire dataset (to save time)
simdvalues <- function(y, ymax, ymin, n, nmax, ci, ntaxa, obsd, currdelta, nsims, lambdavect)  {
   
  # simulate data from sets of thetas corresponding to this delta
  simd <-  rep(NA, nsims)           # vector of simulated values
  
  # initialize thetahiValues to save thetahi values (for testing only)
  thetahiValues <- rep(NA, nsims)
  
  for(sim in 1:nsims)  {     
    thetas <- rep(NA, ntaxa)                      # vector of hypothesized thetas
    
    # choose a random value for the highest theta
    scalefac <- max(ci) - ymax
    lambdaYmax <- lambdavect[which.max(y)]    # lambda for the taxon w/ highest find
    thetahi <- ymax + simulategap(lambdaYmax,nmax)*scalefac    # will this work if there is a tie for ymax?
    thetahi <- max(thetahi, currdelta)        # make sure thetahi is large enough so we can subtract delta
    thetahiValues[sim] <- thetahi             # save thetahi value 

   #  ***** above does not account for the fact that the highest taxon is an order statistic *****
   
    # randomly assign thetahi to one of the taxa
    eligibletaxa <- (1:ntaxa)
    taxonhi <- sample(eligibletaxa, 1)
    thetas[taxonhi] <- thetahi
    
    # set the lowest theta to be currdelta units lower
    thetalo <- thetahi - currdelta

    # randomly assign thetalo to one of the remaining taxa
    remainingtaxa <- (1:ntaxa)[-taxonhi]
    taxonlo <- sample(remainingtaxa, 1)
    thetas[taxonlo] <- thetalo
    
    # set other thetas to intermediate values
    for(i in (1:ntaxa)[-c(taxonlo, taxonhi)]) 
      thetas[i] <- runif(1, thetalo, thetahi)
    
    # simulate dataset from this set of thetas 
    simy <- rep(NA,ntaxa)        # vector of simulated highest finds
    for(i in 1:ntaxa){
      simy[i] <- max(rrefbeta(n[i], lambdavect[i])) * thetas[i]
    }
    # save d from this dataset
    simd[sim] <- max(simy) - min(simy)
  }
  #return(simd)
  return(list(simd, thetahiValues))  # return thetahiValues for testing purposes
}




# Calculate a CI for the range of extinction durations consistent with the data
deltaCI41 <- function(x, conf=.9, PLOT=0, yplotmax=NULL, xplotmax=NULL, 
                      nsims1=100000, nsims2=1000, givenlambdas=NULL, mindiff=5)
{
  # x = dataset; assumes each taxon (column) is sorted high to low and padded with NAs
  # conf = CI confidence level
  # PLOT = make plot?
  # yplotmax, xplotmax = axis limits of plot
  # nsims1 = no. of simulations to run at each value of delta to get sampling dist. for lower endpoint
  # nsims2 = no. of simulations to run at each value of delta to get sampling dist. for upper endpoint
  # givenlambdas = vector of lambdas to pass in (for testing)
  # mindiff = how close the binary search steps need to be to halt
  
    
  # for debugging line-by-line
#  x <- data; conf=.9; PLOT=1; yplotmax=NULL; xplotmax=NULL; nsims1=100000; 
#  nsims2=1000; givenlambdas=NULL; mindiff=5


  
  # ----------------- Initialize and get dataset parameters ----------------- #

  ntaxa <- dim(x)[2]                                               # number of taxa
  n <- rep(NA, ntaxa)                       # vector of sample sizes for each taxon
  for(i in 1:ntaxa) 
    n[i] <- length(na.omit(x[,i]))                # get sample sizes for each taxon
  y <- x[1,]                        # vector of highest fossil finds for each taxon
  ymax <- max(y);                 ymin <- min(y);    
  nmax <- n[which.max(y)]    # sample size of the highest taxon (not the largest sample size)
  # ***** above won’t work right if more than one taxon is tied for the max *****
  obsd <- ymax - ymin                                            # observed d value
  ci <- rep(NA, ntaxa)                  # vector of estimated thetas for each taxon
  ci2 <- rep(NA, ntaxa)                 # vector of estimated thetas for each taxon at higher conf. level
  lambdavect <- rep(NA, ntaxa)          # vector of estimated lambdas for each taxon
  lambdavars <- rep(NA, ntaxa)          # vector of estimated variances for the lambdas of each taxon
  if(is.null(xplotmax))    xplotmax <- obsd*3               # default x-axis limits 
  if(is.null(yplotmax))    yplotmax <- ymax*2               # default y-axis limits 
  thetahiValues <- c()                               # table to save thetahi values
  
  
  # calculate range extensions on each taxon using Adaptive Beta Method (Wang et al 2016)
  # using confidence level 'conf' (e.g., 90%)
  for(i in 1:ntaxa) { 
    temp <- abm38(na.omit(x[,i]), conf=conf)
    ci[i] <- temp[3]                        # upper limit of CI for each taxon
    lambdavect[i] <- temp[4]                # estimated lambda values for each taxon
    lambdavars[i] <- temp[5]                # estimated variances of lambda values for each taxon
  }
  
  
  # set up plot
  if(PLOT)  {
    par(mfrow=c(2,1), mai=c(1,1,.2,.5))
    plotrangechart(x, ci)
    plot(0,0, type="n", xlim=c(0,xplotmax), ylim=c(0,yplotmax), xlab="delta ",
         ylab="simulated d values")    
    abline(h=obsd, lty=3, col=gray(.7))
    abline(v=obsd, lty=3, col=gray(.7))
  }
  
  
  
    
  # ------------- Shrink lambdas towards the mean of all lambdas ------------- #

  # first calculate mean and variance of the lambdas 
  meanLambdas <- mean(lambdavect)
  varLambdas <- var(lambdavect)
  
  # calculate scaling factor to shrink lambdas towards the grand mean
  b <- lambdavars / (lambdavars + rep(varLambdas, ntaxa))
  
  # shrink the lambdas
  lambdaShrunkVect <- b*meanLambdas + (1-b)*lambdavect
  lambdaShrunkVect[lambdaShrunkVect>0] <- 0    # lambdas assumed <= 0; otherwise, shrink to 0
  
  if (!is.null(givenlambdas))      # if lambdas are passed in rather than estimated (for testing)
    lambdaShrunkVect <- givenlambdas
  
  

  
  # ---------- Calculate lower CI endpoint using Fay et al 2015 method ---------- #

  m <- 100
  A <- runif(m)
  B <- runif(m)
  U1 <- U2 <- rep(NA, ntaxa)
  alpha <- (1-conf)/2

  # get the taxa with the highest and lowest CI tops (assumes no ties)
  mintaxon <- which.min(ci)
  maxtaxon <- which.max(ci)
  xmintaxon <- na.omit(x[ , mintaxon])
  xmaxtaxon <- na.omit(x[ , maxtaxon])
  
  # calculate CI endpoints for lower taxon
  L1 <- y[mintaxon]
  for(i in 1:m)
    U1[i] <- abm38(xmintaxon, conf=A[i])[3]
  
  # calculate CI endpoints for higher taxon
  L2 <- y[maxtaxon]
  for(i in 1:m)
    U2[i] <- abm38(xmaxtaxon, conf=B[i])[3]

  # calculate CI endpoints for the difference
  deltamin <- as.numeric( quantile(L2 - U1, alpha) )
  deltamax <- as.numeric( quantile(U2 - L1, 1-alpha) )  # not using this as the actual upper endpoint, just a loop limit


  # ------------- Calculate upper endpoint of confidence interval ------------- #

  low <- obsd           # assumes obsd must be in the CI, so upper endpoint >= obsd
  high <- deltamax      # upper endpoint from Fay method, which is conservative
  currdelta <- obsd
  difference <- high - low
  whichStep <- 0

  # loop starting at currdelta = obsd, iterate until we converge on the upper endpoint
  while(difference > mindiff)  {  

    # simulate d values and determine cutoff for acceptance
    whichStep <- whichStep + 1    
    temp <- simdvalues(y, ymax, ymin, n, nmax, ci, ntaxa, obsd, currdelta, nsims2, lambdaShrunkVect)
    simd <- unlist(temp[1])
    thetahiValues <- c(thetahiValues, unlist(temp[2]))     # ***** just for testing; delete it in final version
    cutoff <- quantile(simd, (1-conf)/2)          # lower-tail 1-sided rejection region
    
    # check if our observed d value is accepted (i.e., consistent w/ the simulated d values)
    if(obsd >= cutoff) {      # if so, search to the right
      if(PLOT)  {
        points(rep(currdelta,nsims2), simd, col=gray(.8), cex=.2)
        lines(rep(currdelta,2), quantile(simd,c(.25,.75)), col="blue", lwd=4) 
        lines(rep(currdelta,2), quantile(simd,c(.05,.95)), col="blue", lwd=2) 
        points(currdelta, cutoff, pch="-", col="blue")
        mtext(whichStep, side=3, line=0, at=currdelta, col="blue", cex=.6)
      }
      low <- currdelta      # make the current delta the new lower bound, keep the old upper bound
    } 
    else {    # if not accepted (observed d value not consistent with simulated values), search to the left
      if(PLOT)  {
        points(rep(currdelta,nsims2), simd, col=gray(.8), cex=.2)
        lines(rep(currdelta,2), quantile(simd,c(.25,.75)), col="steel blue 2", lwd=4) 
        lines(rep(currdelta,2), quantile(simd,c(.05,.95)), col="steel blue 2", lwd=2) 
        points(currdelta, cutoff, pch="-", col="steel blue 2")
        mtext(whichStep, side=3, line=0, at=currdelta, col="steel blue 2", cex=.6)
      }
      high <- currdelta     # make the current delta the new upper bound, keep the old lower bound
    }

    # set the new delta value and update the difference between the lower and upper search bounds
    currdelta <- (low + high)/2
    difference = high - low

  }  # while
  # when we exit this while loop, we have converged to the upper edge of the acceptance region
  # note if mindiff is too small, we may halt before the edge is found
  
  # save upper endpoint
  deltamax <- currdelta
  
  
  
  

  # ---------- adjust endpoints ---------- #
  
  # adjust enpoints of CI to predict abs(max(theta) - min(theta)) instead of max(theta) - min(theta)
  # assumes upper endpoint is already positive
  deltamax <- max(deltamax, abs(deltamin))
  deltamin <- max(deltamin, 0)




  
  # ------------------------------ Plot results ------------------------------ #
  
  # plot CI endpoints on plot
  if(PLOT)  {
     abline(v=deltamin, col="pink")
     abline(v=deltamax, col="steel blue 2")
  }
  
  
  
  # -------- Return CI lower and upper bound, and observed value of d -------- #

  return(c(deltamin,deltamax, obsd))

}


