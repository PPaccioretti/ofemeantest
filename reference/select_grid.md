# Filter grid by minimum observations and minimum number of treatments

Spatially joins points to the grid, counts observations and distinct
treatments per cell, and keeps only cells that satisfy both thresholds.

## Usage

``` r
select_grid(grid, data, x, min_per_cell = 1L, return_points = FALSE)
```

## Arguments

- grid:

  An \`ofe_grid\` (from \`make_grid()\`) or an \`sf\` polygons grid with
  \`CellID\`.

- data:

  An \`sf\` object of points.

- x:

  Character scalar with the name of the treatment column in \`data\`.

- min_per_cell:

  Integer, minimum number of observations per cell (default 1).

- return_points:

  Logical, if \`TRUE\` also returns the joined points restricted to the
  selected cells (component \`points_sel\`). Default \`FALSE\`.

## Value

An \`ofe_grid\` (list) with updated \`grid_sel\` and \`params\`. Adds:

- \`cell_stats\`: a data.frame with per-cell \`n_obs\`, \`n_trt\`, and
  \`treatments\` list.

- \`points_sel\` (optional): the points that fall in the selected cells.
