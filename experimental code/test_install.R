
# List of required packages
req_list=c("foreach", "doParallel")

# Check if devtools is installed
if (!require(devtools)) install.packages("devtools");

# Installing package "linearsdr" from github;
devtools::install_github("HarrisQ/linearsdr");


#### Setting up Parallelization #### 
# Find out how many cores are available (if you don't already know)
print( paste( as.character(detectCores()), "cores detected" ) );
# Create cluster with desired number of cores
cl <- makePSOCKcluster(detectCores()-1)
# Register cluster
doParallel::registerDoParallel(cl)
# Find out how many cores are being used
print( paste( as.character(getDoParWorkers() ), "cores registered" ) )
## [1] 3


# Generate some discrete data


# Test OPCG, RADE, MADE with plots and distance



# Generate some continuous data
# set.seed()


# Test OPCG, RADE, MADE with plots and distance



# End 
remove.packages("linearsdr")