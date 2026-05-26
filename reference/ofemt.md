# OFE permutation analysis

Runs a cell-based OFE analysis with permutation ANOVA and spatial
diagnostics.

## Usage

``` r
ofemt(
  data,
  y,
  x,
  cellsize = 10,
  min_per_cell = 4,
  nmin_cell = NULL,
  n_p = 1000,
  n_s = 200,
  alpha = 0.05,
  shift = c(0, 0),
  p_adjust_method = c("none", "bonferroni", "holm", "BH"),
  crs = NULL,
  keep_components = c("none", "light", "full"),
  grid = NULL,
  angle_deg = 0,
  buffer = 0,
  seed = 7L,
  alpha_bonferroni = NULL
)
```

## Arguments

- data:

  sf points; must contain columns \`y\` and \`x\`.

- y:

  response column (numeric).

- x:

  treatment column (factor/character).

- cellsize, shift, angle_deg, buffer, min_per_cell:

  grid settings (used when \`grid\` is NULL).

- nmin_cell, alpha_bonferroni:

  \*\*Deprecated.\*\* Use \`min_per_cell\` instead of \`nmin_cell\`, and
  \`p_adjust_method = "bonferroni"\` instead of \`alpha_bonferroni\`.

- n_p:

  number of permutations per ANOVA run.

- n_s:

  number of sampling runs.

- alpha:

  significance threshold for letters.

- p_adjust_method:

  p-value adjustment method across pairwise comparisons.

- crs:

  optional target projected CRS if \`data\` is in lon/lat.

- keep_components:

  \`"none"\`, \`"light"\`, or \`"full"\` to embed geometries in output.

- grid:

  optional: an \`ofe_grid\` (e.g., from \[make_ofe_grid()\]). If
  \`NULL\`, the grid is generated internally using \`cellsize\`,
  \`min_per_cell\`, \`shift\`, \`angle_deg\` and \`buffer\`, and an
  informative message is emitted.

- seed:

  integer; controls reproducibility of the permutation sampling (the
  grid itself is deterministic from its construction arguments). Set
  \`seed = NULL\` to let results vary across runs.

## Value

An object of class \`ofemt_result\`.

a list of length 4 with mean test comparisson results

## Details

The methodology involves:

1.  calculation of effective sample size (ESS) given the underlying
    spatial structure.

2.  ANOVA permutation test on a random sample of ESS.

3.  generation of the empirical distribution of p-values from repetition
    of step two. The median of this empirical distribution is regarded
    as the p-value associated with the non-treatment effect hypothesis.

The test can be easily extended to cover scenarios with more than two
treatments by employing pairwise comparisons of OFE across multiple
treatments, with p-values can be adjusted for multiplicity using
Bonferroni correction

## References

A new method to compare treatments in unreplicated on-farm
experimentation. Córdoba M., Paccioretti P., Balzarini M. Under review.

## Examples

``` r
if (FALSE) { # \dontrun{
  my_data <- ofe_f2
  res <- ofemt(my_data, y = "Yield", x = "Treatment",
               cellsize = 10, min_per_cell = 4, alpha = 0.05)
} # }
```
