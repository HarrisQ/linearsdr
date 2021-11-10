
# Linear Sufficient Dimension Reduction

The ‘linearsdr’ package contains popular methods for sufficient
dimension reduction as well as some more recent methods in the
forthcoming paper “Generalized Forward Sufficient Dimension Reduction
for Categorical and Ordinal Responses”.

# Installation

The package can be installed by running:

<!-- ::: .pkgdown-devel -->

``` r
# Install development version from GitHub
devtools::install_github("HarrisQ/linearsdr")
```

<!-- ::: -->

The ‘NAMESPACE - no oxy’ is the original NAMESPACE file prior to when I
started using rOxygen and is kept as a reference. The ‘NAMESPACE - oxy’
is NAMESPACE file generated from rOxygen.

## Current State of Package:

  - Things that work:
      - Forward Linear SDR methods: OPG/OPCG, MAVE/MADE, Tuning OPCG  
      - Inverse Linear SDR methods: SIR, DR, SAVE with options for
        regularization  
  - Minor Things to do:
      - Clean up OPG/OPCG and MAVE/MADE code  
      - (Maybe) Limit function exports to just the ones people will
        use  
      - (Maybe) Clean Regularized OPG/OPCG (RADE) Code
  - Major Things to do:
      - Write Vignettes/Examples  
      - Finish code for Multivariate Continuous Response  
      - Write more Documentation for functions
