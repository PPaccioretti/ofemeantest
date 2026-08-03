# Changelog

## ofemeantest 0.0.900.9000

### Breaking changes

- **[`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
  no longer accepts `nmin_cell` or `alpha_bonferroni`.** Both were kept
  as deprecated aliases; they are now removed. Use `min_per_cell`
  instead of `nmin_cell`, and `p_adjust_method = "bonferroni"` instead
  of `alpha_bonferroni = TRUE`. Since the package has never been on
  CRAN, the disruption should be limited to a rename in existing
  scripts. Note that dropping the two arguments also shifts the
  *positional* order of the remaining ones, so calls that relied on
  argument position (rather than name) need to be revisited.
- **`p_adj` is now the median of the per-run adjusted p-values.**
  Previously the median p-value was computed first and adjusted
  afterwards. Now
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) is
  applied within each sampling run, across that run’s pairwise
  comparisons, and the reported `p_adj` is the median of the resulting
  distribution. This keeps `p_value` and `p_adj` on the same footing —
  both medians over the `n_s` runs — and makes the median line in
  [`plot_pvalue_hist()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_pvalue_hist.md)
  agree with the table. Numeric results change when
  `p_adjust_method != "none"`.
- `keep_components = "light"` now stores the complete `ofe_grid` object
  (full grid, selection, per-cell stats and parameters) rather than only
  the selected cells, so the result can be plotted directly.

### Bug fixes

- Treatment labels are no longer rewritten internally. Spaces used to be
  replaced by dots (`"Bortrac + Zintrac"` became `"Bortrac.+.Zintrac"`),
  which leaked into `Means comparison`, `Cells per treatment` and the
  `Comparison` column, and could produce `NA` letters when a label
  collided with the separator used by
  [`multcompView::multcompLetters()`](https://lselzer.github.io/multcompView/reference/multcompLetters.html).
  Labels are now mapped to internal placeholders for the whole
  computation and restored on output, so any label — spaces, `+`, `-`,
  parentheses, accents — is reported verbatim. The same rewriting is
  gone from
  [`select_grid()`](https://ppaccioretti.github.io/ofemeantest/reference/select_grid.md),
  where it could collapse two genuinely different treatments into one
  and make a mixed cell look single-treatment.
- The compact letter display now follows the ordering of the means:
  treatments are sorted by decreasing median response before the letters
  are assigned, so `"a"` always marks the highest-yielding group.
- The upper bound passed to
  [`spdep::dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.html)
  is exclusive, so the one cell whose nearest-neighbour distance
  *defined* `dmax` was being dropped from its own neighbour graph. The
  bound is now nudged to include it. This silently degraded `Rho`,
  Moran’s I and the ESS on every analysis; on the bundled `ofe_f2`
  example the corrected values are ESS 97 / Rho 0.652 / Moran’s I 0.548
  (previously 99 / 0.646 / 0.529).
- [`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
  no longer flips `spdep`’s global `ZeroPolicyOption` as a side effect
  when an empty neighbour set is encountered; the policy is passed
  explicitly to the calls that need it.
- Response and treatment column names are quoted before being used in a
  formula, so columns whose names contain spaces or other special
  characters (`"Rinde kg/ha (seco)"`) no longer fail. Point coordinates
  are stored under dot-prefixed names, so user columns literally named
  `X` or `Y` are no longer shadowed.
- Sampling a treatment that contributes a single cell no longer falls
  into [`sample()`](https://rdrr.io/r/base/sample.html)’s length-one
  trap (`sample(5, 1)` sampling from `1:5`).

### Reporting and plotting

- `spdep`’s *“neighbour object has N sub-graphs”* warnings are muffled.
  They are expected for OFE trial shapes and carry no consequence for
  the diagnostics. What is now reported instead is the condition that
  actually matters: an error when no cell has any neighbour (the
  distance threshold is too small for the selection), and a warning
  naming how many cells were left isolated.
- [`plot_grid_selection()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_grid_selection.md)
  accepts an `ofemt_result` directly, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods are
  provided for both `ofe_grid` and `ofemt_result`. On a result produced
  with `keep_components = "full"` a bare `plot(res)` draws the full
  grid, the selected cells and the observations, with no further
  arguments — the view needed to judge the effect of `cellsize`,
  `shift`, `angle_deg` and `buffer`. Selected cells are shaded
  translucently so the points underneath stay visible, and the plot
  title repeats the grid parameters so successive calls can be compared.
- Both plotting functions gained an `engine` argument, `"ggplot2"`
  (default) or `"base"`. The ggplot2 engine places the legend *outside*
  the plotting panel, so it can never sit on top of the data and the
  result no longer depends on the size of the graphics device — the base
  engine has to fit the legend inside the panel, which on a dense grid
  or a small window left it unreadable. `ggplot2` stays in `Suggests`:
  when it is not installed both functions fall back to `"base"` with a
  message.
- [`plot_grid_selection()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_grid_selection.md)
  gained `point_size` (and `legend_pos`, which accepts base-style corner
  keywords under either engine). Yield-monitor data runs to tens of
  thousands of observations, where the default dot size still reads as a
  solid mass; lowering it brings the cell boundaries back.
- Unselected grid cells are now filled a very light grey instead of left
  transparent, so they read as cells rather than as page background and
  their legend key is visible.
- [`plot_pvalue_hist()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_pvalue_hist.md)
  now draws a solid line at the median of the distribution in addition
  to the dashed line at `alpha`, and defaults to the adjusted p-values
  whenever the analysis used a `p_adjust_method` (override with
  `which = "raw"` or `which = "adjusted"`).
- [`print()`](https://rdrr.io/r/base/print.html) reports the total
  number of grid cells alongside the number selected, states which
  adjustment produced `p_adj`, and notes that the means table is sorted
  by decreasing mean.
- [`?ofemt`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
  now documents every component of the returned object — including
  `perm_runs`, the per-run empirical p-value distribution the method is
  built on — and spells out what each `keep_components` level stores and
  costs.

### Other changes

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
