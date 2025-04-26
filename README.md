
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ovbsa

<!-- badges: start -->
<!-- badges: end -->

The goal of `ovbsa` is to conduct sensitivity analysis of omitted
variable bias in linear econometrics models.

## Installation

You can install the development version of ovbsa from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("dbasu-umass/ovbsa")
```

## Example 1

This is a basic example which shows you how to find the bias-adjusted
estimate, the bias-adjusted standard error and the bias-adjusted
confidence interval in a linear regression model. The user needs to
choose a benchmark covariate, values of the sensitivity parameters `kD`
and `kY` and the significance level `alpha` for testing the null
hypothesis that the treatment effect is zero.

Let us first load the relevant libraries.

``` r
library(ovbsa)
library(sensemakr)
#> See details in:
#> Carlos Cinelli and Chad Hazlett (2020). Making Sense of Sensitivity: Extending Omitted Variable Bias. Journal of the Royal Statistical Society, Series B (Statistical Methodology).
```

In this example we will use use the data set `darfur` from the package
`sensemakr`, which studies the effect of exposure to violence on
attitudes towards peace. Here we choose `kD=3`, `kY=3` and `alpha=0.05`

``` r
## basic example code
analysis1 <- ovbsa::bsal(
  kd=3,ky=3,alpha=0.05,data=darfur,
  outcome = "peacefactor",
  treatment = "directlyharmed",
  bnch_reg = "female",
  other_reg = c("village","age","farmer_dar","herder_dar",
                "pastvoted","hhsize_darfur")
)
```

Now, let us see the results.

``` r
(analysis1)
#>                        Case 1       Case 2       Case 3
#> r2yd.x            0.021873093  0.021873093  0.021873093
#> r2dz.x            0.008040002  0.027492860  0.027492860
#> r2yz.dx           0.781414454  0.380969988  0.374050471
#> estimate          0.097315819  0.097315819  0.097315819
#> adjusted_estimate 0.045525612  0.029779889  0.030396023
#> adjusted_se       0.010924114  0.018566585  0.018670065
#> adjusted_lower_CI 0.024081595 -0.006666285 -0.006253282
#> adjusted_upper_CI 0.066969628  0.066226063  0.067045329
```

## Example 2

Continuing with the previous example, we will now compute the
probability that taking account of omitted variable bias will overturn
the conclusion of the study. Here we choose `alpha=0.05`.

First, total R-squared based analysis:

``` r
# total r2-based analysis
res1 <- ovbsa::saltr2(
  alpha=0.05,data = darfur, outcome = "peacefactor",
  treatment = "directlyharmed", bnch_reg = "female",
  other_reg = c("village","age","farmer_dar","herder_dar",
                "pastvoted","hhsize_darfur"),
  N = 500, k_kd=1, k_ky=1
)
# probability
(res1$frac_prob_wt)
#> [1] 0.2582741
```

Second: partial R-squared based analysis without conditioning on the
treatment variable:

``` r
res2 <- ovbsa::salpr2ncd(
  alpha=0.05,data = darfur, outcome = "peacefactor",
  treatment = "directlyharmed", bnch_reg = "female",
  other_reg = c("village","age","farmer_dar","herder_dar",
                "pastvoted","hhsize_darfur"),
  N = 500, k_kd=1, k_ky=1
)
# probability
(res2$frac_prob_wt)
#> [1] 0.3118923
```

Finally, partial R-squared based analysis with conditioning on the
treatment variable:

``` r
res3 <- ovbsa::salpr2cd(
  alpha=0.05,data = darfur, outcome = "peacefactor",
  treatment = "directlyharmed", bnch_reg = "female",
  other_reg = c("village","age","farmer_dar","herder_dar",
                "pastvoted","hhsize_darfur"),
  N = 500, k_kd=1, k_ky=1
)
# probability
(res3$frac_prob_wt)
#> [1] 0.3384226
```
