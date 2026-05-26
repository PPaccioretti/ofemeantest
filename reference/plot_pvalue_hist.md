# Plot histogram(s) of permutation p-values per comparison

Produces one or several histograms showing the distribution of
permutation-derived \*p\*-values for each comparison contained in an
\[\`ofemt_result\`\] object. A dashed red vertical line marks the
significance threshold (\`alpha\`) used in the analysis.

## Usage

``` r
plot_pvalue_hist(results)
```

## Arguments

- results:

  An object of class \`ofemt_result\`, typically obtained from
  \[ofemt()\]. Must contain \`perm_runs\` (data frame with
  \`Comparison\` and \`p_value\`) and \`params\$alpha\`.

## Value

A \`ggplot\` object if \*\*ggplot2\*\* is installed; otherwise the
function falls back to base R \`hist()\` and returns \`NULL\` invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment", cellsize = 9)
  plot_pvalue_hist(res)
} # }
```
