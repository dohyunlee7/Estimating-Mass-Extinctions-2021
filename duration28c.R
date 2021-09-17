# ---------------------------------------------------------------- #
#
# Confidence interval for extinction duration
# Steve Wang, Phil Everson, Brendan McVeigh, Aaron Zimmerman, Heidi Wong
# from Wang et al 2012, Paleobiology
# Confidence intervals for the duration of a mass extinction
# R code, version 28c
#
# The main function is deltaCI, at the end. 
#
# ---------------------------------------------------------------- #




# Strauss and Sadler uniform recovery CI
rangeext <- function(x, conflevel=.9)  {
  x <- na.omit(x)
  xmax <- max(x);      n <- length(x)
  return( xmax*(1-conflevel)^(-1/n) )
}





# Draw range chart
plotrangechart <- function(data, ci)  {
  ord <- order(data[1,])
  data <- data[,ord]           # sort data by highest find
  ci <- ci[ord]
  ntaxa <- dim(data)[2]
  maxfinds <- dim(data)[1]
  ymax <- max(data[1,])
  plot( NULL, pch="", xlab="", ylab="", xlim=c(0,ntaxa+1), 
        ylim=c(-.3,ymax*2.1), bty="L", xaxt="n" )
  for(taxon in 1:ntaxa)  {
    # draw taxon lines
    lines( c(taxon,taxon), c(0,data[1,taxon]), lwd=.9, col="darkgray" )
    # plot fossil horizons
    points(rep(taxon,maxfinds), data[,taxon], pch=16, cex=.65)
    # plot confidence interval tops
    points(taxon, ci[taxon], pch="^", cex=.65, col="red")
  }
}




 


simdvalues <- function(y, ymax, ymin, n, nmax, nmin, ci, ntaxa, obsd, currdelta, nsims)  {
# simulate a dataset from thetas having true delta equal to ‘delta’, 
# and calculate the d value from the dataset.
# differs from simdata() in that only the highest finds are simulated (to save time)

  # Simulate data from sets of thetas corresponding to this delta
  simd <-  rep(NA,nsims)                          # vector of simulated values
  for(sim in 1:nsims)  {     
    thetas <- rep(NA, ntaxa)                      # vector of hypothesized thetas

    # fix thetahi and set thetalo to be currdelta units lower
    # choose a random value for the highest theta
    scalefac <- max(ci) - ymax
    thetahi <- ymax + rbeta(1,1,nmax)*scalefac 
    thetahi <- max(thetahi, currdelta)
   
    # randomly assign thetahi to one of the taxa
    eligibletaxa <- (1:ntaxa)
    taxonhi <- sample( eligibletaxa, 1 )
    thetas[taxonhi] <- thetahi

    # set the lowest theta
    thetalo <- thetahi - currdelta
    # randomly assign thetalo to one of the remaining taxa
    remainingtaxa <- (1:ntaxa)[-taxonhi]
    taxonlo <- sample( remainingtaxa, 1 )
    thetas[taxonlo] <- thetalo

    # set other thetas to intermediate values
    for(i in (1:ntaxa)[-c(taxonlo, taxonhi)]) 
      thetas[i] <- runif(1, thetalo, thetahi)

    # simulate dataset from this set of thetas 
    simy <- rep(NA,ntaxa)        # vector of simulated highest finds
    for(i in 1:ntaxa)
      simy[i] <- max(runif(n[i], 0, thetas[i]))

    # save d from this dataset
    simd[sim] <- max(simy) - min(simy)
  }

  return(simd)
}



 



deltaCI28c <- function(x, conflevel=.9, stepsize=4, PLOT=0, yplotmax=0, xplotmax=0, nsims=300)
{
  # Calculate a CI for the range of extinction durations consistent with the data
  # x = dataset; assumes each taxon is sorted high to low and padded with NAs
  # conflevel = CI confidence level
  # stepsize = resolution to loop over
  # nsims = no. of simulations to run at each value of delta to get sampling dist.


  # Initialize and get dataset parameters
  ntaxa <- dim(x)[2]
  n <- rep(NA, ntaxa)                       # vector of sample sizes for each taxon
  for(i in 1:ntaxa)
    n[i] <- length(na.omit(x[,i]))                # get sample sizes for each taxon
  y <- x[1,]                        # vector of highest fossil finds for each taxon
  ymax <- max(y);                 ymin <- min(y);    
  nmax <- n[which.max(y)];        nmin <- n[which.min(y)]
          # This won’t work right if more than one taxon is tied for the max or min
  obsd <- ymax - ymin                                            # observed d value
  if(xplotmax==0)    xplotmax <- obsd+6*stepsize            # default x-axis limits 
  if(yplotmax==0)    yplotmax <- ymax*1.7                   # default y-axis limits 


  # Calculate range extensions on each taxon
  ci <- rep(NA, ntaxa) 
  for(i in 1:ntaxa)  
    ci[i] <-rangeext(x[,i], conflevel=conflevel)

  # Set up plot
  if(PLOT)  {
    par(mfrow=c(2,1), mai=c(1,1,.2,.5))
    plotrangechart(x, ci)
    plot(0,0, type="n", xlim=c(-stepsize,xplotmax), ylim=c(0,yplotmax), xlab="delta ",
         ylab="simulated d values")    
    abline(h=obsd, col=gray(.6))
  }


  # Calculate lower endpoint of confidence interval
  # Loop starting at delta = zero until we accept
  accept <- 0                                                    # do we accept here?
  currdelta <- -stepsize                                                 # initialize
  while(accept==0)  {

    # simulate data from sets of thetas corresponding to this delta
    currdelta <- currdelta + stepsize                    # begin with currdelta = 0
    simd <- simdvalues(y, ymax, ymin, n, nmax, nmin, ci, ntaxa, obsd, currdelta, nsims)

    # check if our observed d value is consistent w/ the simulated d values
    cutoff <- quantile(simd, 1-(1-conflevel)/2)          # 1-sided rejection region
    if(obsd <= cutoff)  
      accept <- 1                   # observed value consistent with simulated values
    if(PLOT)  {
      points(rep(currdelta,nsims), simd, col=gray(.8), cex=.2)
      lines(rep(currdelta,2), quantile(simd,c(.25,.75)), col="red", lwd=4) 
      lines(rep(currdelta,2), quantile(simd,c(.05,.95)), col="red", lwd=2) 
    }
  } 
  #  # save accepted value, which is the lower endpoint
  deltamin <- currdelta - stepsize/2
  deltamin <- max(deltamin, 0)                                  # can’t go below zero
  # as we exit this loop, accept = 1 and currdelta = deltamin + stepsize/2
  # i.e., we are in the acceptance region. It is assumed that currdelta < obsd.
  
  

  # Calculate upper endpoint of confidence interval
  # Loop going up until we reject; then everything higher must reject   
  while(accept)  {

    # simulate data from sets of thetas corresponding to this delta
    currdelta <- currdelta + stepsize    
    simd <- simdvalues(y, ymax, ymin, n, nmax, nmin, ci, ntaxa, obsd, currdelta, nsims)

    # check if our observed d value is consistent w/ the simulated d values
    cutoff <- quantile(simd, (1-conflevel)/2)                # 1-sided rejection region
    if(obsd < cutoff)  
      accept <- 0             # observed value not consistent with simulated values
    if(PLOT)  {
      color <- ifelse(currdelta < obsd, "purple", "blue")
      points(rep(currdelta,nsims), simd, col=gray(.8), cex=.2)
      lines(rep(currdelta,2), quantile(simd,c(.25,.75)), col=color, lwd=4) 
      lines(rep(currdelta,2), quantile(simd,c(.05,.95)), col=color, lwd=2) 
    }
  } 
  # if we exit the loop, then currdelta was rejected; save the previous currdelta
  deltamax <- currdelta - stepsize/2    

  # Return CI lower and upper bound, and observed value of d
  if(PLOT)    points(c(deltamin,deltamax), c(0,0), pch=6, col=gray(.7))
  return(c(deltamin,deltamax, obsd))
}

