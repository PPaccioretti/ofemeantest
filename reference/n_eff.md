# Effective sample size under spatial autocorrelation

Heuristic function to adjust sample size \`n\` by a spatial
autocorrelation parameter \`rho\`. Based on: Griffith, D.A., &
Peres-Neto, P.R. (2006). Spatial Modeling in Ecology: The Flexibility of
Eigenfunction Spatial Analyses. Ecology, 87(10), 2603–2613.

## Usage

``` r
n_eff(n, rho)
```

## Arguments

- n:

  integer, nominal sample size.

- rho:

  numeric in \[0, 1\], spatial autocorrelation intensity.

## Value

numeric, effective sample size.
