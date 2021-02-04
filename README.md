# linearsdr
 Linear Sufficient Dimension Reduction methods

Contains code for my forthcoming paper on "Generalized Forward Regression for Sufficient Dimension Reduction with Multinomial Response". 

The code can be installed as a package via the install_github function in devtools:
`devtools::install_github("HarrisQ/linearsdr");`

The 'NAMESPACE - no oxy' is the original NAMESPACE file prior to when I started using rOxygen and is kept as a reference.
The 'NAMESPACE - oxy' is NAMESPACE file generated from rOxygen.

The 'Experimental Code' folder is where I have stored code for aspects of Linear SDR that I have been working on sparingly and is not quite ready for use. 

Things that work:
---
* OPG/OPCG, MAVE/MADE  
* SIR, DR, SAVE  
* OPCG Tuning 

Minor Things to do:
---
* Clean up OPP/OPCG and MAVE/MADE code  
* Limit function exports to just the ones peopel will use  

Major Things to do:
---
* Build Test Script  
* Finish/Clean Code for OPCG and Multivariate Continuous Response  
* Finish/Clean Regularization (RADE) Code  
* Write more Documentation  




