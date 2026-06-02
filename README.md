
<!-- README.md is generated from README.Rmd. Please edit that file -->

# melsvmp

<!-- badges: start -->

<!-- badges: end -->

The package “melsvmp” is an implementation of a Variational Message
Passing algorithm for fitting Mixed-Effect Location-Scale model, which
is defined as:

``` math
\begin{aligned}
y_{ij} &= x_{ij}^\top \beta + \nu_i + \varepsilon_{ij} \\
\nu_i &\sim \mathcal{N}(0, \exp{(U_i^\top \alpha)}) \\
\varepsilon_{ij} &\sim \mathcal{N}(0, \exp{(W_{ij}^\top \tau + \omega_i)}) \\
\omega_i &\sim \mathcal{N}(0, \sigma_\omega^2)
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
a MELS model. We use an EMA dataset “posmood” to illustrate the usage,
you can import and view the description of this dataset by

``` r
library(melsvmp)

data("posmood")
```

The following code shows how to fit a MELS model with our vmp algorithm
on the posmood dataset. The standard errors and confidence intervals are
estimated by the sandwich estimators proposed by Westling and McCormick
(2019).

``` r
posmood_vmp <- mels_vmp(y = "posmood", 
                            beta_formula = ~ other_bs + other_ws + genderf + t1 + t2 +
                              t3 + t4 + w1 + w2 + w3 + w4 + w5 + w6 + tirbor + frustr, 
                            alpha_formula = ~ other_bs + genderf + age15, 
                            tau_formula = ~ other_bs + other_ws + genderf + age15 + tirbor + frustr, 
                            id = "id", 
                            data = posmood)
#> CAVI converges at iteration  18
summary(posmood_vmp)
#> ## VMP for MELS (with Robust Sandwich Errors) ##
#> --------------------------------------------------------
#> --- Mean Model Parameters (beta) ---
#>             Estimate Robust SE CI.Lower CI.Upper  z value p-value
#> (Intercept)   8.1434    0.2664   7.6212   8.6655  30.5674  0.0000
#> other_bs      0.4616    0.3262  -0.1777   1.1010   1.4152  0.1570
#> other_ws      0.1911    0.0288   0.1347   0.2475   6.6455  0.0000
#> genderf      -0.0291    0.0944  -0.2142   0.1560  -0.3083  0.7578
#> t1            0.2472    0.0381   0.1726   0.3217   6.4944  0.0000
#> t2            0.3550    0.0427   0.2714   0.4387   8.3220  0.0000
#> t3            0.3841    0.0410   0.3038   0.4644   9.3738  0.0000
#> t4            0.2998    0.0590   0.1842   0.4155   5.0817  0.0000
#> w1            0.1403    0.0460   0.0501   0.2304   3.0481  0.0023
#> w2            0.1686    0.0519   0.0668   0.2704   3.2470  0.0012
#> w3            0.1278    0.0527   0.0246   0.2311   2.4271  0.0152
#> w4            0.0821    0.0476  -0.0112   0.1754   1.7253  0.0845
#> w5            0.0502    0.0429  -0.0340   0.1343   1.1685  0.2426
#> w6           -0.0219    0.0406  -0.1014   0.0576  -0.5395  0.5896
#> tirbor       -0.1457    0.0081  -0.1615  -0.1298 -18.0511  0.0000
#> frustr       -0.3262    0.0093  -0.3443  -0.3080 -35.1695  0.0000
#> 
#> --- Between-Subject Variance Parameters (alpha) ---
#>             Estimate Robust SE CI.Lower CI.Upper z value p-value
#> (Intercept)  -0.1878    0.3498  -0.8733   0.4978 -0.5369  0.5914
#> other_bs      0.3324    0.4775  -0.6036   1.2683  0.6960  0.4864
#> genderf       0.0926    0.1456  -0.1928   0.3780  0.6360  0.5248
#> age15        -0.1040    0.0652  -0.2318   0.0238 -1.5950  0.1107
#> 
#> --- Within-Subject Variance Parameters (tau) ---
#>             Estimate Robust SE CI.Lower CI.Upper z value p-value
#> (Intercept)   0.1931    0.1709  -0.1419   0.5281  1.1298  0.2586
#> other_bs     -0.3058    0.2074  -0.7124   0.1007 -1.4744  0.1404
#> other_ws     -0.0921    0.0347  -0.1600  -0.0241 -2.6563  0.0079
#> genderf       0.1680    0.0575   0.0553   0.2808  2.9205  0.0035
#> age15        -0.0850    0.0259  -0.1358  -0.0343 -3.2845  0.0010
#> tirbor       -0.0042    0.0079  -0.0197   0.0113 -0.5346  0.5929
#> frustr        0.1300    0.0078   0.1147   0.1452 16.7232  0.0000
#> 
#> --- Random Effect Standard Deviation (omega) ---
#>               Estimate
#> omega_std_dev   0.5882
```

For performing percentile bootstrap, you can use the following code, we
support both parallel and sequential computing. This can give a more
robust confidence interval.

``` r
# posmood_vmp_boot <- bootstrap_mels_vmp(posmood_vmp, B = 1000, cores = 10)
posmood_vmp_boot <- bootstrap_mels_vmp(posmood_vmp, B = 1000, parallel = FALSE)
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
summary(posmood_vmp_boot)
#> ## Bootstrap Summary for MELS Model ##
#> --------------------------------------
#> Successful replicates: 1000
#> Total runtime: 1.15 hours
#> 
#> --- Mean Model Parameters (beta) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)   8.1434  0.2663   7.5945   8.6468
#> other_bs      0.4616  0.3285  -0.1582   1.1190
#> other_ws      0.1911  0.0299   0.1329   0.2496
#> genderf      -0.0291  0.0914  -0.2073   0.1527
#> t1            0.2472  0.0381   0.1774   0.3264
#> t2            0.3550  0.0433   0.2763   0.4448
#> t3            0.3841  0.0408   0.3097   0.4659
#> t4            0.2998  0.0594   0.1922   0.4211
#> w1            0.1403  0.0478   0.0511   0.2360
#> w2            0.1686  0.0516   0.0727   0.2695
#> w3            0.1278  0.0534   0.0260   0.2419
#> w4            0.0821  0.0485  -0.0110   0.1756
#> w5            0.0502  0.0426  -0.0309   0.1333
#> w6           -0.0219  0.0404  -0.0969   0.0568
#> tirbor       -0.1457  0.0082  -0.1606  -0.1289
#> frustr       -0.3262  0.0095  -0.3450  -0.3083
#> 
#> --- Between-Subject Variance Parameters (alpha) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)  -0.1878  0.3791  -1.0010   0.5178
#> other_bs      0.3324  0.5205  -0.6398   1.4259
#> genderf       0.0926  0.1452  -0.1985   0.3754
#> age15        -0.1040  0.0688  -0.2364   0.0294
#> 
#> --- Within-Subject Variance Parameters (tau) ---
#>             Estimate Boot.SE CI.Lower CI.Upper
#> (Intercept)   0.1931  0.1748  -0.1342   0.5364
#> other_bs     -0.3058  0.2120  -0.6913   0.1331
#> other_ws     -0.0921  0.0363  -0.1621  -0.0219
#> genderf       0.1680  0.0566   0.0617   0.2859
#> age15        -0.0850  0.0264  -0.1339  -0.0308
#> tirbor       -0.0042  0.0082  -0.0208   0.0126
#> frustr        0.1300  0.0078   0.1154   0.1452
#> 
#> --- Random Effect Standard Deviation (omega) ---
#>               Estimate Boot.SE CI.Lower CI.Upper
#> omega_std_dev   0.5882   0.023    0.539   0.6283
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
