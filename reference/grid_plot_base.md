# Draw the grid selection with base graphics

Draw the grid selection with base graphics

## Usage

``` r
grid_plot_base(
  grid_all,
  grid_sel,
  data,
  main,
  legend,
  legend_pos,
  point_size = NULL,
  ...
)
```

## Arguments

- grid_all, grid_sel:

  \`sf\` polygons; \`grid_sel\` may be \`NULL\`.

- data:

  Optional \`sf\` points to overlay. Rarely needed: the observations are
  taken from the object itself — \`points_sel\` for an \`ofe_grid\`
  built with \`make_ofe_grid(return_points = TRUE)\`, \`points_joined\`
  for an \`ofemt_result\` run with \`keep_components = "full"\`. Pass
  \`data\` when the object carries no points, or to override the stored
  ones.

- main:

  Plot title. Defaults to a one-line summary of the grid parameters,
  which is what makes successive calls comparable.

- legend:

  Logical; draw the legend. Default \`TRUE\`.

- legend_pos:

  Where to place the legend. Defaults to the engine's own sensible
  choice: \`"right"\` (outside the panel) for \`"ggplot2"\`, and
  \`"topleft"\` for \`"base"\`. Base-style keywords are translated for
  the ggplot2 engine, so \`"bottomright"\` works with either; \`"none"\`
  hides it.

- point_size:

  Size of the observation dots. Defaults to \`0.15\` for the ggplot2
  engine and \`0.35\` (as \`cex\`) for the base engine. Yield-monitor
  data runs to tens of thousands of points, where the default can still
  read as a solid mass — lower it to see the cell boundaries underneath.

- ...:

  Further arguments passed to the underlying \[plot()\] call for the
  full-grid layer. Base engine only; ignored by the ggplot2 engine.

## Value

Invisibly \`NULL\`.
