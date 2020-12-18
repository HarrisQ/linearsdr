
# Check if devtools is installed
if (!require(devtools)) install.packages("devtools");

# Installing package "linearsdr" from github;
devtools::install_github("HarrisQ/linearsdr");

# Generate some discrete data


# Test OPCG, RADE, MADE with plots and distance


# Generate some continuous data


# Test OPCG, RADE, MADE with plots and distance



# End 
remove.packages("linearsdr")