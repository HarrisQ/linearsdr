# Linear Sufficient Dimension Reduction (linearsdr)
 Linear Sufficient Dimension Reduction methods

Contains code for my forthcoming paper on "Generalized Forward Regression for Sufficient Dimension Reduction with Multinomial Response". 

The code can be installed as a package via the install_github function in devtools:
`devtools::install_github("HarrisQ/linearsdr");`

The 'NAMESPACE - no oxy' is the original NAMESPACE file prior to when I started using rOxygen and is kept as a reference.
The 'NAMESPACE - oxy' is NAMESPACE file generated from rOxygen.

The 'Experimental Code' folder is where I have stored code for aspects of Linear SDR that I have been working on sparingly and is not quite ready for use. 

Current State of Package:
----
* Things that work:
    + Forward Linear SDR methods: OPG/OPCG, MAVE/MADE, Tuning OPCG  
    + Inverse Linear SDR methods: SIR, DR, SAVE with options for regularization  
* Minor Things to do:
    + Clean up OPG/OPCG and MAVE/MADE code  
    + (Maybe) Limit function exports to just the ones people will use  
    + (Maybe) Clean Regularized OPG/OPCG (RADE) Code  

* Major Things to do:
    + Write Vignettes/Examples   
    + Finish code for Multivariate Continuous Response   
    + Write more Documentation for functions  




