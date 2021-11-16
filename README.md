
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Linear Sufficient Dimension Reduction

<!-- # pkgdown <img src="man/figures/logo.png" align="right" alt="" width="120" /> -->

The ‘linearsdr’ package contains popular methods for sufficient
dimension reduction as well as some more recent methods in the
forthcoming paper “Generalized Forward Sufficient Dimension Reduction
for Categorical and Ordinal Responses”.

## Installation

The package can be installed by running:

<!-- ::: .pkgdown-devel -->

``` r
# Install development version from GitHub
devtools::install_github("HarrisQ/linearsdr")
```

<!-- ::: -->

## Current State of Package:

  - Current features of the package:
      - Inverse Linear SDR methods: SIR, DR, SAVE with options for
        regularization
      - Forward Linear SDR methods: OPG, MAVE, OPCG, MADE  
      - Some vignettes demonstrating the use of methods in package
  - Immediate tasks to complete:
      - Finish code for Multivariate Continuous Response for inverse
        methods
      - Write more Documentation for functions
  - Not-so pressing tasks to complete:
      - Clean up OPG/OPCG and MAVE/MADE code  
      - Limit function exports to just the ones people will use  
      - Speed up MAVE’s Code
        <!-- + Clean Regularized OPG/OPCG (RADE) Code   -->
