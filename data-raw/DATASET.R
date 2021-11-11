## code to prepare `airfoil_datta` dataset goes here

airfoil_datta = read.table("data-raw/airfoil_self_noise.dat", header=F, skip=0)


usethis::use_data(airfoil_datta, overwrite = TRUE)
