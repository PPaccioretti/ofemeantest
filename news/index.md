# Changelog

## ofemeantest 0.0.900.9000

- [`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
  exposes a `seed` argument that controls reproducibility of the
  permutation sampling. Default is `seed = 7L`; pass `NULL` to let
  results vary across runs.
- [`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
  accepts a pre-built grid via the `grid` argument. When `grid = NULL`
  the function builds one internally with
  [`make_ofe_grid()`](https://ppaccioretti.github.io/ofemeantest/reference/make_ofe_grid.md);
  the two paths are interchangeable when called with matching arguments.
- [`make_ofe_grid()`](https://ppaccioretti.github.io/ofemeantest/reference/make_ofe_grid.md)
  exported as the canonical way to construct and inspect the analysis
  grid before running the test.
- [`plot_grid_selection()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_grid_selection.md)
  and
  [`plot_pvalue_hist()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_pvalue_hist.md)
  added for visual checks of the grid selection and the permutation
  p-value distribution.
  [`plot_pvalue_hist()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_pvalue_hist.md)
  has a base R fallback when `ggplot2` is not installed.
- [`print()`](https://rdrr.io/r/base/print.html) method for
  `ofemt_result` objects provides a concise console summary.
- Two new vignettes:
  [`vignette("getting-started")`](https://ppaccioretti.github.io/ofemeantest/articles/getting-started.md)
  (user manual) and
  [`vignette("methodology")`](https://ppaccioretti.github.io/ofemeantest/articles/methodology.md)
  (technical description of the protocol).
- `pkgdown` website at <https://ppaccioretti.github.io/ofemeantest/>.
- Dependencies reduced: dropped `dplyr`, `magrittr` and `tidyr` from
  `Imports`. `ggplot2` moved to `Suggests`.
- [`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
  rewritten on top of the new grid pipeline; cells are filtered once in
  [`select_grid()`](https://ppaccioretti.github.io/ofemeantest/reference/select_grid.md),
  not twice.
- Moran’s I is now extracted by name (`$estimate["Moran I statistic"]`)
  instead of by position, with a dependency-contract test guarding
  against future renames in `spdep`.
- Test suite added (`testthat` 3rd edition) covering grid construction,
  end-to-end
  [`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
  on `ofe_f2`, reproducibility from `seed`, error paths, and contracts
  with upstream dependencies.
- Continuous integration on GitHub Actions: `R CMD check` on
  Ubuntu/macOS/Windows times R release+devel+oldrel; pkgdown site
  deployed on push to `main`.
- `nmin_cell` is deprecated; use `min_per_cell` instead.
- `alpha_bonferroni` is deprecated; use `p_adjust_method = "bonferroni"`
  instead.
