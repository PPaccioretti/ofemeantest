# Print method for \`ofemt_result\` objects

Displays a concise summary of an \[\`ofemt_result\`\] object, including
general information about the experimental grid, descriptive statistics,
and permutation-based test results.

## Usage

``` r
# S3 method for class 'ofemt_result'
print(x, ...)
```

## Arguments

- x:

  An object of class \[\`ofemt_result\`\], typically returned by
  \[ofemt()\]. The object must contain, at minimum, a \`"General
  information"\` data frame. Optionally, it may include components such
  as \`"Means comparison"\` and \`"ANOVA permutation test"\`.

- ...:

  Additional arguments passed to or from other methods (ignored).

## Value

Invisibly returns the input object \`x\`, unchanged.

## Details

This print method provides a human-readable summary in the console:

- \*\*General information:\*\* Cell size, number of selected cells, and
  per-cell observation counts (minimum, median, maximum).

- \*\*Spatial statistics:\*\* Sample size (\`n\`), effective sample size
  (\`ESS\`), spatial correlation (\`Rho\`), and Moran’s \*I\*.

- \*\*Means comparison:\*\* Shows the first 10 rows of the table of
  comparisons between treatment means.

- \*\*ANOVA permutation test:\*\* Displays pairwise tests summarizing
  median \*p\*-values and corrected \*p\*-values across permutation
  runs.

If any of these components are missing, they are simply skipped in the
printed output.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- ofemt(my_data, y = "Yield", x = "Treatment")
  print(res)
} # }
```
