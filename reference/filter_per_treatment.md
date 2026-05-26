# Filter OFE grid cells containing a single treatment

Performs a spatial join between point data (\`sf\`) and an OFE grid,
identifies cells containing exactly one unique treatment, and returns an
updated \`ofe_grid\` object with only those eligible cells retained in
its \`grid_sel\` component.

## Usage

``` r
filter_per_treatment(data, grid, x)
```

## Arguments

- data:

  An \`sf\` object of points containing the treatment column \`x\`.

- grid:

  An \`ofe_grid\` object (as returned by \[make_grid()\] or similar).
  Its \`grid_sel\` component will be used and replaced with the filtered
  version. Must be polygonal and share CRS with \`data\`.

- x:

  A character string of length one giving the name of the treatment
  column in \`data\`.

## Value

An updated \`ofe_grid\` object whose \`grid_sel\` component contains
only polygons corresponding to single-treatment cells.

## Details

The function:

1.  Validates the input types and column name.

2.  Replaces spaces in treatment labels with dots and stores them in a
    temporary column \`".trt"\`.

3.  Ensures that the grid has a \`"CellID"\` column (creates one if
    missing).

4.  Joins points to grid cells via \`sf::st_intersects()\`.

5.  Aggregates per cell and keeps only cells containing a single unique
    treatment (\`unique_trt == 1\`).

Points outside the grid are dropped (\`left = FALSE\` in the join).

## Examples

``` r
if (FALSE) { # \dontrun{
  filtered_grid <- filter_per_treatment(data = pts_sf, grid = my_ofe_grid, x = "Treatment")
  plot(filtered_grid$grid_sel["CellID"])
} # }
```
