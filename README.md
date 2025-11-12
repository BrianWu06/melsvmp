
<!-- README.md is generated from README.Rmd. Please edit that file -->

# melsvmp

<!-- badges: start -->

<!-- badges: end -->

The goal of melsvmp is to …

## Installation

You can install the development version of melsvmp like so:

``` r
install.packages("devtools")
#> Installing package into 'C:/Users/user/AppData/Local/Temp/RtmpSaIK2o/temp_libpath29747183190c'
#> (as 'lib' is unspecified)
#> package 'devtools' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\user\AppData\Local\Temp\RtmpuIZy3c\downloaded_packages
devtools::install_github("BrianWu06/melsvmp")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo BrianWu06/melsvmp@HEAD
#> rbibutils  (2.3   -> 2.4  ) [CRAN]
#> reformulas (0.4.1 -> 0.4.2) [CRAN]
#> Installing 2 packages: rbibutils, reformulas
#> Installing packages into 'C:/Users/user/AppData/Local/Temp/RtmpSaIK2o/temp_libpath29747183190c'
#> (as 'lib' is unspecified)
#> package 'rbibutils' successfully unpacked and MD5 sums checked
#> package 'reformulas' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\user\AppData\Local\Temp\RtmpuIZy3c\downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>       ✔  checking for file 'C:\Users\user\AppData\Local\Temp\RtmpuIZy3c\remotes460c4c587dbf\BrianWu06-melsvmp-49efb55/DESCRIPTION'
#>       ─  preparing 'melsvmp':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>       ─  building 'melsvmp_0.0.0.9000.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/user/AppData/Local/Temp/RtmpSaIK2o/temp_libpath29747183190c'
#> (as 'lib' is unspecified)
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
#>   Total Runtime: 0.15 secs
```

``` r
riesby_vmp_boot <- bootstrap_mels_vmp(riesby_vmp, B = 1000, cores = 10)
#> Starting PARALLEL bootstrap with 1000 replicates on 10 cores...
#> Warning: All bootstrap replicates failed.
#> Warning: All bootstrap replicates failed.
#> Warning: All bootstrap replicates failed.
#> Warning: All bootstrap replicates failed.
summary(riesby_vmp_boot)
#> ## Bootstrap Summary for MELS Model ##
#> --------------------------------------
#> Successful replicates: 0
#> Total runtime: 2.98 secs
#> 
#> --- Mean Model Parameters (beta) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)  22.2573      NA       NA       NA
#> week         -2.2673      NA       NA       NA
#> endog         1.8679      NA       NA       NA
#> endweek      -0.0139      NA       NA       NA
#> 
#> --- Between-Subject Variance Parameters (alpha) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)   2.2777      NA       NA       NA
#> endog         0.4937      NA       NA       NA
#> 
#> --- Within-Subject Variance Parameters (tau) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)   2.1381      NA       NA       NA
#> week          0.1841      NA       NA       NA
#> endog         0.2928      NA       NA       NA
#> 
#> --- Random Effect Standard Deviation (omega) ---
#>               Estimate Boot.SE CI.Lower CI.Upper
#> omega_std_dev   0.6101      NA       NA       NA
#> --------------------------------------
```
