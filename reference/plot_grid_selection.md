# Plot the grid, the selected cells and the observations

Draws the three layers that matter when tuning a grid, on one set of
axes: the \*\*full grid\*\* in light grey, the \*\*selected cells\*\*
(those that passed the single-treatment and \`min_per_cell\` filters)
shaded and outlined in black, and the \*\*observations\*\* as points.
Seeing the points against the cell boundaries is the quickest way to
judge the effect of \`cellsize\`, \`shift\`, \`angle_deg\` and
\`buffer\`.

## Usage

``` r
plot_grid_selection(
  x,
  data = NULL,
  points = TRUE,
  main = NULL,
  legend = TRUE,
  legend_pos = "topleft",
  ...
)

# S3 method for class 'ofe_grid'
plot(x, ...)

# S3 method for class 'ofemt_result'
plot(x, ...)
```

## Arguments

- x:

  Either an \`ofe_grid\` (from \[make_ofe_grid()\]) or an
  \`ofemt_result\` (from \[ofemt()\]) that was run with
  \`keep_components = "light"\` or \`"full"\`. With an \`ofemt_result\`
  every layer is taken from the object itself, so no extra arguments are
  needed; the points are only available under \`keep_components =
  "full"\`.

- data:

  Optional \`sf\` points to overlay. Rarely needed: the observations are
  taken from the object itself — \`points_sel\` for an \`ofe_grid\`
  built with \`make_ofe_grid(return_points = TRUE)\`, \`points_joined\`
  for an \`ofemt_result\` run with \`keep_components = "full"\`. Pass
  \`data\` when the object carries no points, or to override the stored
  ones.

- points:

  Logical; set to \`FALSE\` to skip the point layer. When \`TRUE\` (the
  default) and no observations are available, a message explains how to
  obtain them rather than silently drawing a grid without points.

- main:

  Plot title. Defaults to a one-line summary of the grid parameters,
  which is what makes successive calls comparable.

- legend:

  Logical; draw the legend. Default \`TRUE\`.

- legend_pos:

  Where to place the legend, passed to \[graphics::legend()\] (e.g.
  \`"topleft"\`, \`"bottomright"\`, or \`"top"\`). Default
  \`"topleft"\`.

- ...:

  Further arguments passed to the underlying \[plot()\] call for the
  full-grid layer.

## Value

Invisibly returns \`NULL\`. Called for its side effect (a base R plot).

## See also

\[make_ofe_grid()\], \[ofemt()\]

## Examples

``` r
if (FALSE) { # \dontrun{
  g <- make_ofe_grid(ofe_f2, x = "Treatment", cellsize = 9, min_per_cell = 4)
  plot_grid_selection(g, data = ofe_f2)
  plot(g, data = ofe_f2)          # same thing

  res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment", cellsize = 9,
               keep_components = "full")
  plot(res)                        # grid + selection + points, no extra args
} # }
```
