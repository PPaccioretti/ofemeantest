# Plot histogram(s) of permutation p-values per comparison

Produces one histogram per pairwise comparison showing the empirical
distribution of permutation \*p\*-values across the \`n_s\` sampling
runs of an \[\`ofemt_result\`\]. Two reference lines are drawn on each
panel: a solid line at the \*\*median\*\* of the distribution — the
value \[ofemt()\] reports for that comparison — and a dashed line at the
significance threshold \*\*alpha\*\*.

## Usage

``` r
plot_pvalue_hist(results, which = c("auto", "adjusted", "raw"), bins = 30)
```

## Arguments

- results:

  An object of class \`ofemt_result\`, typically obtained from
  \[ofemt()\]. Must contain \`perm_runs\` and \`params\$alpha\`.

- which:

  Which p-values to plot. \`"auto"\` (the default) uses the
  multiplicity-adjusted values whenever the analysis was run with
  \`p_adjust_method != "none"\`, and the raw values otherwise.
  \`"adjusted"\` and \`"raw"\` force one or the other. The adjustment is
  applied within each run before the histogram is built, so the median
  line coincides with the \`p_adj\` reported in the \`ANOVA permutation
  test\` table.

- bins:

  Number of histogram bins. Default 30.

## Value

A \`ggplot\` object if \*\*ggplot2\*\* is installed; otherwise the
function falls back to base R \`hist()\` and returns \`NULL\` invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment", cellsize = 9,
               p_adjust_method = "bonferroni")
  plot_pvalue_hist(res)              # adjusted p-values
  plot_pvalue_hist(res, which = "raw")
} # }
```
