---
title: "Linear SDR for multivariate Y"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Linear SDR for multivariate Y}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.0.5
library(linearsdr)
data('energy_datta', package="linearsdr")
# summary(energy_datta)
# dim(energy_datta)

X = as.matrix(energy_datta[,1:8])
Y1 = as.numeric(energy_datta$Y1)
Y2 = as.numeric(energy_datta$Y2)

pairs(X)
```

![](C:/Users/Harri/AppData/Local/Temp/Rtmp6nB301/preview-5eb41576cd4.dir/sdr_energy_files/figure-html/setup-1.png)<!-- -->


```r

b_hat_sir = sir(x=X, y=Y1, nslices = 10, d=2, ytype = "continuous" )$beta


linearsdr:::ggplot_fsdr(Y1, t((X)%*%b_hat_sir[,1:2]), y_on_axis=F,
                        ytype="continuous",
                        h_lab='SIR 1', v_lab='SIR 2',
                        main_lab= paste0('SIR'), size=3)
```

![](C:/Users/Harri/AppData/Local/Temp/Rtmp6nB301/preview-5eb41576cd4.dir/sdr_energy_files/figure-html/unnamed-chunk-2-1.png)<!-- -->



```r

b_hat_dr = dr(x=X, y=Y1, nslices = 10, d=2, ytype = "continuous" )$beta

linearsdr:::ggplot_fsdr(Y1, t((X)%*%b_hat_dr[,1]), y_on_axis=T,
                        ytype="continuous",
                        h_lab='DR 1', v_lab='DR 2',
                        main_lab= paste0('DR'), size=3)
```

![](C:/Users/Harri/AppData/Local/Temp/Rtmp6nB301/preview-5eb41576cd4.dir/sdr_energy_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

Forward methods using OPG

```r
# Parallelization not run because of computational time. 

library("doParallel")
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
library("foreach")
print( paste( as.character(detectCores()), "cores detected" ) );
#> [1] "16 cores detected"
# Create cluster with desired number of cores
cl <- makePSOCKcluster(detectCores()-1)
# Register cluster
doParallel::registerDoParallel(cl)
# Find out how many cores are being used
print( paste( as.character(getDoParWorkers() ), "cores registered" ) )
#> [1] "15 cores registered"
# stopCluster(cl)

X_std=(sapply(1:dim(X)[2], FUN= function(j)
  center_cpp(X[,j], NULL) ) )%*%matpower_cpp(cov((X)) , -1/2);

b_hat_opg = opcg(x=X_std, y=Y1, bw = .75, d=2, ytype = "continuous", 
                 method= "cg", parallelize = T )


linearsdr:::ggplot_fsdr(Y1, t((X_std)%*%b_hat_opg[,1]), y_on_axis=T,
                        ytype="continuous",
                        h_lab='OPG 1', v_lab='OPG 2',
                        main_lab= paste0('OPG'), size=3)
```

![](C:/Users/Harri/AppData/Local/Temp/Rtmp6nB301/preview-5eb41576cd4.dir/sdr_energy_files/figure-html/unnamed-chunk-4-1.png)<!-- -->
