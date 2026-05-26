# Plot grid selection (base R)

Quick visual check of an \`ofe_grid\` object: full grid in grey,
selected cells (those that passed the single-treatment and
\`min_per_cell\` filters) outlined in black. Optionally overlays the
point data.

## Usage

``` r
plot_grid_selection(ofe_grid, data = NULL)
```

## Arguments

- ofe_grid:

  An object of class \`ofe_grid\` (typically returned by
  \[make_ofe_grid()\]).

- data:

  Optional \`sf\` points overlaid on the grid.

## Value

Invisibly returns \`NULL\`. Called for its side effect (a base R plot).
