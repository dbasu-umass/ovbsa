
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ovbsa

<!-- badges: start -->

<!-- badges: end -->

The goal of `ovbsa` (omitted variable bias sensitivity analysis) is to
conduct sensitivity analysis of reported results in linear econometrics
models to the presence of omitted variable bias. The key function in
this package computes a bias-adjusted confidence interval by
implementing the method presented in Basu (2026).

Basu (2026)’s method builds on and extends Cinelli and Hazlett (2020).
Consider a linear econometric model where you believe there is an
unobserved confounder. Choose a regressor that is included in the model
as a benchmark covariate. Let $k_D$ and $k_Y$ denote sensitivity
parameters.

Here $k_D$ captures the strength of association of the unobserved
confounder with the treatment variable (i.e., regressor of interest)
relative to the strength of association of the benchmark covariate with
the treatment variable. Similarly, $k_Y$ captures the strength of
association of the unobserved confounder with the outcome variable
relative to the strength of association of the benchmark covariate with
the outcome variable.

Cinelli and Hazlett (2020) had presented a method to compute omitted
variable bias-adjusted confidence intervals for specific values of $k_D$
and $k_Y$. Basu (2026) asks us to consider $k_D$ and $k_Y$ as random
variables, imposes reasonable prior distributions on these random
variables and then computes bias-adjusted confidence intervals by
integrating out $k_D$ and $k_Y$.

In this package, I use independent truncated exponential distributions
as prior distribution of $k_D$ and $k_Y$. For more details see Basu
(2026).

## Installation

You can install the package `ovbsa` from CRAN with:

``` r
# uncomment this line
# install.packages("ovbsa")
```

You can install the development version of ovbsa from
[GitHub](https://github.com/) with:

``` r
# uncomment these lines
# install.packages("pak")
# pak::pak("dbasu-umass/ovbsa")
```

## Main function

The main function in this package is:

- `baci`: this function computes the bias-adjusted confidence interval
  (baci) in a linear regression model when $k_D$ and $k_Y$ follow
  truncated exponential distributions; the function returns a list
  `results` containing the unadjusted and bias-adjusted confidence
  intervals; it also returns a list `undstats` containing the estimate,
  the unadjusted stamdard error, the max of $k_D$ and the max of $k_Y$;
  it also returns a `ggplot2` plot object visualizing the support of the
  joint distribution of $(k_D,k_Y)$.

Let us first load the relevant libraries and then work through an
example.

``` r
library(ovbsa)
library(sensemakr)
#> See details in:
#> Carlos Cinelli and Chad Hazlett (2020). Making Sense of Sensitivity: Extending Omitted Variable Bias. Journal of the Royal Statistical Society, Series B (Statistical Methodology).
library(ggplot2)
```

In the examples we will use use the data set `darfur` from the package
`sensemakr`, which studies the effect of exposure to violence on
attitudes towards peace.

## Bias-adjusted confidence interval

This is a basic example which shows you how to find the bias-adjusted
estimate, the bias-adjusted standard error and the bias-adjusted
confidence interval in a linear regression model. To use this function
the user needs to choose a benchmark covariate, the significance level
`alpha` for testing the null hypothesis that the treatment effect is
zero and `N` to create a grid of size $N \times N$ over the support of
the distribution of $(k_D,k_Y)$.

Here we choose the benchmark covariate as `female`, `alpha=0.05` and
`N=1000`. In the first step, we estimate the model with ordinary least
squares using the `lm` function.

``` r
## fit model
fit <- lm(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
             pastvoted + hhsize_darfur + female + village,
             data = sensemakr::darfur)
```

In the second step, we call `baci` to conduct sensitivity analysis.

``` r
# conduct sensitivity analysis
res1 <- ovbsa::baci(fit = fit, treatment = "directlyharmed",
                benchmark = "female", N = 1000, alpha = 5/100)
#> Extracting regression quantities : 0.8 sec elapsed
#> Computing values on the grid : 0.36 sec elapsed
#> Prior distribution: truncated exponential : 1.05 sec elapsed
#> Results available now!
```

Let us see the results, starting with the underlying statistics:

``` r
res1$undstats
#>                 Values
#> estimate    0.09731582
#> std error   0.02325654
#> max(kD)   109.11924031
#> max(kY)     8.17145822
```

Let us now see the unadjusted and bias-adjusted confidence intervals.

``` r
res1$results
#>                       Lower     Upper
#> Unadjusted CI    0.05166327 0.1429684
#> Bias-adjusted CI 0.02596471 0.1234509
```

Finally, let us see a picture of the support of the joint prior
distribution of $k_D$ and $k_Y$.

``` r
print(res1$support_kdky_plot)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" alt="" width="100%" />

## References

- Basu, D. (2026). How likely is it that omitted variable bias will
  overturn your results? SSRN Working Paper. Available here:
  <doi:10.2139/ssrn.4704246>

- Cinelli, C. and Hazlett, C. (2020). Making Sense of Sensitivity:
  Extending Omitted Variable Bias. *Journal of the Royal Statistical
  Society Series B: Statistical Methodology,* 82(1):39–67
