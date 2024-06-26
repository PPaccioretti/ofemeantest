
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ofemeantest

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of ofemeantest from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pkg_install("PPaccioretti/ofemeantest")
```

## Example

The example dataset comprised a field trial of corn (Zea mays L.), where
the impact of broadcast phosphorus fertilization was evaluated. The
treatment strip received superphosphate at rates ranging from 250 to 500
$kg\ ha^{–1}$ \[0-21-0\], while the control strip received no P
fertilizer (0 $kg\ ha^{–1}$). Each strip was 2.2 ha in area, totaling
100 ha. Raw grain yield monitor data for each each single strip were
cleaned following the protocol proposed by Vega et al. (2019).

``` r
library(ofemeantest)
data("ofe_f2")
```

The `ofemt` function first creates a regular artificial grid over the
area, including the strips to be analyzed. The grid spacing (`cellsize`)
was set to 9 m. This grid overlaid the point data from the yield map,
assigning each georeferenced data point to a corresponding grid cell.
Each grid cell contained a minimum of four data points
(`nmin_cell = 4`). Empty cells or cells with less than five data points
are removed. The median of the yield data is calculated for each cell. A
one-way ANOVA is fitted to the cell medians to compare treatments. From
the residuals of the fitted model, the magnitude of spatial
autocorrelation is estimated using the approximate profile likelihood
estimator of spatial dependence (Rho) and the Moran index (MI). An
effective sample size (ESS) is then calculated. Second, `ofemt` randomly
samples grid cells without replacement using the ESS. The yield data
from both strips (treatment and control) are used to perform a
permutational ANOVA with 2000 iterations per test (`n_p`). This
procedure provides a p-value for the comparison of treatment means. The
permutational ANOVA is repeated on 200 random samples of grid cells
(`n_s`) using the same sampling procedure. The significance level for
comparing means (`alpha`), is left at their default value (0.05). The
resulting permutation p-values are collected to generate an empirical
p-value distribution. Finally, the median p-value is calculated to
compare the means.

``` r
ofemt(data = ofe_f2,
      y = "Yield_tn",
      x = "Treatment",
      cellsize = 9,
      nmin_cell = 4,
      n_p = 2000,
      n_s = 200,
      alpha = 0.05,
      shift = 0,
      alpha_bonferroni = FALSE,
      crs = NULL)
#> $`General information`
#>   Cellsize Min.Obs.Cell Max.Obs.Cell Median.Obs.Cell   n ESS      Rho        MI
#> 1        9            4            7               5 554  97 0.652285 0.5483085
#> 
#> $`Cells per treatment`
#> 
#>    Control Fertilized 
#>        273        281 
#> 
#> $`ANOVA permutation test`
#>               Comparison p-value
#> 1 Control vs. Fertilized   0.008
#> 
#> $`Means comparison`
#>    Treatment Yield_tn_mean letters
#> 2 Fertilized      5.280402       b
#> 1    Control      4.760095      a
```

The results show that out of the total number of cells selected (554),
the ESS was 97, with a magnitude of autocorrelation (Rho) of 0.65. The
fertilization treatment showed significant differences compared to the
control (p \< 0.05).

The average yield of the fertilized strip was 5.28 $t\ ha^{-1}$, while
that of the control was 4.76 $t\ ha^{-1}$.

## References

Córdoba M., Paccioretti P., Balzarini M. 2024. A new method to compare
treatments in unreplicated on-farm experimentation. Precis. Agric.
(under review).

Vega A., Córdoba M., Balzarini M. 2019. Protocol for automating error
removal from yield maps. Precis. Agric. 20: 1030–1044.
<https://doi.org/10.1007/s11119-018-09632-8>
