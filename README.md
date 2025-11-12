
<!-- README.md is generated from README.Rmd. Please edit that file -->

# melsvmp

<!-- badges: start -->

<!-- badges: end -->

The goal of melsvmp is to …

## Installation

You can install the development version of melsvmp like so:

``` r
# install.packages("devtools")
# devtools::install_github("BrianWu06/melsvmp")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(melsvmp)

data(riesby)

riesby_vmp <- mels_vmp(y = "hamd", 
                       beta_formula = ~ week + endog + endweek, 
                       alpha_formula = ~ endog, 
                       tau_formula = ~ week + endog, 
                       id = "id",
                       data = riesby)
#> CAVI converges at iteration  55
summary(riesby_vmp)
#> ## VMP for MELS (with Robust Sandwich Errors) ##
#> --------------------------------------------------------
#> --- Mean Model Parameters (beta) ---
#>             Estimate Robust SE CI.Lower CI.Upper z value p-value
#> (Intercept)  22.2573    0.6428  20.9974  23.5173 34.6241  0.0000
#> week         -2.2673    0.2855  -2.8269  -1.7077 -7.9410  0.0000
#> endog         1.8679    1.0036  -0.0992   3.8350  1.8612  0.0627
#> endweek      -0.0139    0.4097  -0.8169   0.7891 -0.0339  0.9729
#> 
#> --- Between-Subject Variance Parameters (alpha) ---
#>             Estimate Robust SE CI.Lower CI.Upper z value p-value
#> (Intercept)   2.2777    0.3218   1.6469   2.9084  7.0775   0.000
#> endog         0.4937    0.3759  -0.2430   1.2304  1.3136   0.189
#> 
#> --- Within-Subject Variance Parameters (tau) ---
#>             Estimate Robust SE CI.Lower CI.Upper z value p-value
#> (Intercept)   2.1381    0.2046   1.7372   2.5391 10.4523  0.0000
#> week          0.1841    0.0549   0.0765   0.2918  3.3513  0.0008
#> endog         0.2928    0.2223  -0.1429   0.7286  1.3171  0.1878
#> 
#> --- Random Effect Standard Deviation (omega) ---
#>               Estimate Approx. SE CI.Lower CI.Upper
#> omega_std_dev   0.6101     0.0544   0.5036   0.7166
#> -------------------------------------------------------
#> Convergence Details:
#>   Algorithm converged in 55 iterations.
#>   Total Runtime: 0.13 secs
```

``` r
# riesby_vmp_boot <- bootstrap_mels_vmp(riesby_vmp, B = 1000, cores = 10)
riesby_vmp_boot <- bootstrap_mels_vmp(riesby_vmp, B = 1000, parallel = FALSE)
#> Starting SEQUENTIAL bootstrap with 1000 replicates...
#>   Running replicate: 100 of 1000
#>   Running replicate: 200 of 1000
#>   Running replicate: 300 of 1000
#>   Running replicate: 400 of 1000
#>   Running replicate: 500 of 1000
#>   Running replicate: 600 of 1000
#>   Running replicate: 700 of 1000
#>   Running replicate: 800 of 1000
#>   Running replicate: 900 of 1000
#>   Running replicate: 1000 of 1000
summary(riesby_vmp_boot)
#> ## Bootstrap Summary for MELS Model ##
#> --------------------------------------
#> Successful replicates: 1000
#> Total runtime: 2.09 mins
#> 
#> --- Mean Model Parameters (beta) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)  22.2573  0.6859  21.0114  23.6263
#> week         -2.2673  0.3122  -2.9267  -1.7150
#> endog         1.8679  1.0300  -0.1190   3.8602
#> endweek      -0.0139  0.4231  -0.8161   0.7724
#> 
#> --- Between-Subject Variance Parameters (alpha) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)   2.2777  0.4671   1.1963   2.7931
#> endog         0.4937  0.5242  -0.2940   1.6136
#> 
#> --- Within-Subject Variance Parameters (tau) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)   2.1381  0.2439   1.6174   2.5498
#> week          0.1841  0.0594   0.0759   0.3012
#> endog         0.2928  0.2545  -0.1872   0.8016
#> 
#> --- Random Effect Standard Deviation (omega) ---
#>               Estimate Boot.SE CI.Lower CI.Upper
#> omega_std_dev   0.6101  0.1372   0.3461   0.8824
#> --------------------------------------
```
