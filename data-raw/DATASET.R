## code to prepare datasets goes here
set.seed(2021)

#### NASA Airfoil Dataset ----
airfoil_datta = read.table("data-raw/airfoil_self_noise.dat", header=F, skip=0)

sample_set = sample(1:dim(airfoil_datta)[1], 500)
airfoil_datta = airfoil_datta[sample_set,]

usethis::use_data(airfoil_datta, overwrite = TRUE)

#### Energy Efficiency Dataset ----
# Bivariate Y

energy_datta = readxl::read_xlsx("data-raw/ENB2012_data.xlsx")

usethis::use_data(energy_datta, overwrite = TRUE)

#### Superconductor Dataset ----

supercond_datta_entire = read.csv2("data-raw/superconductor.csv", header = T, sep=",")
dim(supercond_datta_entire)


sample_set = sample(1:dim(supercond_datta_entire)[1], 500)
supercond_datta = supercond_datta_entire[sample_set,]
supercond_datta = sapply(1:dim(supercond_datta)[2], 
                         function(k) as.numeric(supercond_datta[,k]))

usethis::use_data(supercond_datta, overwrite = TRUE)

#### Communities and Crime Dataset ----

crime_datta_entire=as.data.frame(read.csv('data-raw/communities.data', header = F, sep = ','))
# dim(crime_datta_entire)
# Drop first 5; drop NA; 102-118, 122-125, 127
crime_datta = crime_datta_entire[,-c(1:5, 31,102:118, 122:125, 127)]

sample_set = sample(1:dim(crime_datta)[1], 500)
crime_datta = crime_datta[sample_set,]

crime_datta = sapply(1:dim(crime_datta)[2], 
                     function(k) as.numeric(crime_datta[,k]))

usethis::use_data(crime_datta, overwrite = TRUE)


#### Pen digit ----

set.seed(2020)
pendigit_train_all=as.data.frame(read.csv('data-raw/pendigits_tra.txt', header = F, sep = ','))
pendigit_test_all=as.data.frame(read.csv('data-raw/pendigits_tes.txt', header = F, sep = ','))

train_set=sample(1:dim(pendigit_train_all)[1],1000) #1:dim(uspsdigit_train)[1];#
test_set=sample(1:dim(pendigit_test_all)[1],1000) #1:dim(uspsdigit_train)[1];#

pendigits_datta = rbind(pendigit_train_all[train_set,],
                        pendigit_test_all[test_set,] )

usethis::use_data(pendigits_datta, overwrite = TRUE)

