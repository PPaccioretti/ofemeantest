# Package index

## Main analysis

The entry point: runs the cell-based permutation analysis end-to-end.

- [`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
  : OFE permutation analysis

## Grid construction

Build, inspect and customise the spatial grid before passing it to
[`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md).

- [`make_ofe_grid()`](https://ppaccioretti.github.io/ofemeantest/reference/make_ofe_grid.md)
  : Create an OFE grid and select valid cells
- [`plot_grid_selection()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_grid_selection.md)
  : Plot grid selection (base R)

## Diagnostics and plotting

- [`plot_pvalue_hist()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_pvalue_hist.md)
  : Plot histogram(s) of permutation p-values per comparison
- [`print(`*`<ofemt_result>`*`)`](https://ppaccioretti.github.io/ofemeantest/reference/print.ofemt_result.md)
  : Print method for \`ofemt_result\` objects

## Bundled data

- [`ofe_f2`](https://ppaccioretti.github.io/ofemeantest/reference/ofe_f2.md)
  : On-Farm single strip treatment trial (Field F2)
- [`ofe_p`](https://ppaccioretti.github.io/ofemeantest/reference/ofe_p.md)
  : On-Farm single strip treatment trial (companion dataset)
