
<!-- README.md is generated from README.Rmd. Please edit that file -->

# melsvmp

<!-- badges: start -->

<!-- badges: end -->

The package “melsvmp” is an implementation of a Variational Message
Passing algorithm for fitting Mixed-Effect Location-Scale model, which
is defined as:

``` math
\begin{aligned}
y_{ij}&= x_{ij}^\top \beta + \nu_i + \varepsilon_{ij} \\
\nu_i&\sim \mathcal{N}(0, \exp{(U_i^\top \alpha)}) \\
\varepsilon_{ij}&\sim \mathcal{N}(0, \exp{(W_{ij}^\top \tau + \omega_i)}) \\
\omega_i&\sim \mathcal{N}(0, \sigma_\omega^2)
\end{aligned}
```

For more information, please refer to Hedeker and Nordgren (2013). Note
that our model formulation is a restricted to only allow for non
time-varying covariates for $\alpha$ and without the effect of $\nu_i$
on $\omega_i$.

Some reviews of Variational Inference (VI) and Variational Message
Passing (VMP) algorithms can be found in Blei, Kucukelbir, and McAuliffe
(2017) and Ormerod and Wand (2010). Briefly speaking, VI introduces
surrogate distributions $q(\theta)$ to approximate the posterior
distribution $p(\theta \mid x)$ by minimizing the KL-divergence of these
two distributions. While the KL-divergence is often intractable, we will
instead optimize the Evidence Lower Bound (ELBO):

``` math
\mathrm{ELBO} = \mathbf{E}_q[\log p(\theta, x)] - \mathbf{E}_q[\log q(\theta)]
```

A common approach to make this approximation computationally efficient
is to use the mean-field assumption, which assumes
$q(\Theta) = q(\theta_1) \cdots q(\theta_m)$ for the set of all latent
variables $\Theta = \{ \theta_1, \dots, \theta_m \}$. Under this
assumption, we can get clean update equations for conjugate variables
(Blei, Kucukelbir, and McAuliffe (2017)) and non-conjugate ones using
Laplace approximation (Wand (2014)).

## Installation

You can install the development version of melsvmp like so:

``` r
# install.packages("devtools")
# devtools::install_github("BrianWu06/melsvmp")
```

## Example

This is a basic example showing how to use the function mels_vmp to fit
a MELS model by an example. We use a simple longitudinal dataset
“riesby” to illustrate the usage, you can import and view the
description of this dataset by

``` r
library(melsvmp)

data(riesby)

?riesby
```

The following code shows how to fit a MELS model with our vmp algorithm
on the riesby dataset. The standard errors and confidence intervals are
estimated by the sandwich estimators proposed by Westling and McCormick
(2019).

``` r
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
#>   Total Runtime: 0.4 secs
```

For performing percentile bootstrap, you can use the following code, we
support both parallel and sequential computing. This can give a more
robust confidence interval.

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
#> Total runtime: 3.47 mins
#> 
#> --- Mean Model Parameters (beta) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)  22.2573  0.6638  20.9902  23.6405
#> week         -2.2673  0.3164  -2.9014  -1.6663
#> endog         1.8679  1.0288  -0.1943   3.8257
#> endweek      -0.0139  0.4268  -0.8583   0.8545
#> 
#> --- Between-Subject Variance Parameters (alpha) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)   2.2777  0.5373   1.2321   2.7799
#> endog         0.4937  0.5863  -0.2390   1.5850
#> 
#> --- Within-Subject Variance Parameters (tau) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)   2.1381  0.2425   1.5929   2.5428
#> week          0.1841  0.0612   0.0643   0.3087
#> endog         0.2928  0.2534  -0.1637   0.7974
#> 
#> --- Random Effect Standard Deviation (omega) ---
#>               Estimate Boot.SE CI.Lower CI.Upper
#> omega_std_dev   0.6101  0.1377   0.3441   0.8749
#> --------------------------------------
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-blei2017variational" class="csl-entry">

Blei, David M, Alp Kucukelbir, and Jon D McAuliffe. 2017. “Variational
Inference: A Review for Statisticians.” *Journal of the American
Statistical Association* 112 (518): 859–77.

</div>

<div id="ref-hedeker2013mixregls" class="csl-entry">

Hedeker, Donald, and Rachel Nordgren. 2013. “MIXREGLS: A Program for
Mixed-Effects Location Scale Analysis.” *Journal of Statistical
Software* 52: 1–38.

</div>

<div id="ref-ormerod2010explaining" class="csl-entry">

Ormerod, John T, and Matt P Wand. 2010. “Explaining Variational
Approximations.” *The American Statistician* 64 (2): 140–53.

</div>

<div id="ref-wand2014fully" class="csl-entry">

Wand, Matt P. 2014. “Fully Simplified Multivariate Normal Updates in
Non-Conjugate Variational Message Passing.” *Journal of Machine Learning
Research*.

</div>

<div id="ref-westling2019beyond" class="csl-entry">

Westling, T, and TH McCormick. 2019. “Beyond Prediction: A Framework for
Inference with Variational Approximations in Mixture Models.” *Journal
of Computational and Graphical Statistics* 28 (4): 778–89.

</div>

</div>
