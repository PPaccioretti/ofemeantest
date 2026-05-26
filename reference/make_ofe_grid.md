# Create an OFE grid and select valid cells

Convenience wrapper that builds a grid with \[all_cells_grid()\] and
filters cells using \[select_grid()\]. It is the main entry point for
generating grids for on-farm experiment (OFE) analysis.

## Usage

``` r
make_ofe_grid(
  data,
  x,
  cellsize,
  min_per_cell = 1L,
  angle_deg = 0,
  buffer = 0,
  shift = c(0, 0),
  return_points = FALSE
)
```

## Arguments

- data:

  An \`sf\` object of points. See \[all_cells_grid()\].

- x:

  Character scalar with the treatment column name in \`data\`. See
  \[select_grid()\].

- cellsize:

  Grid cell size. See \[all_cells_grid()\].

- min_per_cell:

  Minimum observations per cell. See \[select_grid()\].

- angle_deg:

  Grid rotation in degrees. See \[all_cells_grid()\].

- buffer:

  Buffer applied to the data bounding box before gridding. See
  \[all_cells_grid()\].

- shift:

  Offset applied to the grid origin. See \[all_cells_grid()\].

- return_points:

  Logical. See \[select_grid()\].

## Value

An \`ofe_grid\` object returned by \[select_grid()\], containing the
full grid, the selected cells, and associated metadata.

## See also

\[all_cells_grid()\], \[select_grid()\]
