# Translate a legend position to ggplot2's vocabulary

Lets the same \`legend_pos\` value work with either engine: base's
corner keywords collapse to the nearest ggplot2 side, and anything
ggplot2 already understands (including a numeric \`c(x, y)\`) passes
through untouched.

## Usage

``` r
legend_pos_gg(pos, legend = TRUE)
```

## Arguments

- pos:

  \`legend_pos\` as supplied by the user; \`NULL\` for the default.

- legend:

  Logical; \`FALSE\` forces \`"none"\`.

## Value

A value suitable for \`ggplot2::theme(legend.position = )\`.
