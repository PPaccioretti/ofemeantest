# Build OFE grid and (initially) select all cells

Build OFE grid and (initially) select all cells

## Usage

``` r
all_cells_grid(data, cellsize, angle_deg = 0, buffer = 0, shift = c(0, 0))
```

## Arguments

- data:

  An \`sf\` object of points with geometry.

- cellsize:

  numeric length-1 or length-2, grid cell size.

- angle_deg:

  numeric, grid rotation in degrees.

- buffer:

  numeric, buffer added around data bbox before gridding (same units as
  CRS).

- shift:

  numeric length-2, offset added to grid origin (xmin, ymin).

## Value

An object of class \`ofe_grid\` (a list) with:

- \`grid_all\`: \`sf\` polygons of the full grid (with \`CellID\`).

- \`grid_sel\`: initially identical to \`grid_all\` (no filtering yet).

- \`params\`: a named list of grid parameters (cellsize, shift, angle,
  buffer, CRS).
