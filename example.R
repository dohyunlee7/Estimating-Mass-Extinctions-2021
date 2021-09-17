# read in a dataset and run the deltaCI function
# you may need to change the file path names 

# read in the code
source("~/Desktop/duration41.4.R")

# read in dataset
data <- matrix(scan("~/Desktop/data11.txt"), 12, 2, byrow=T)  

# run deltaCI function (takes about a minute on my laptop)
deltaCI41(data, PLOT=0)
