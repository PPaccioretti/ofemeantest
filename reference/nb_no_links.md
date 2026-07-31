# Regions with no neighbours in an \`nb\` object

In \`spdep\`, a region with an empty neighbour set is stored as the
integer \`0L\`. This helper returns the positions of those regions.

## Usage

``` r
nb_no_links(nb)
```

## Arguments

- nb:

  An \`nb\` object (e.g. from \`spdep::dnearneigh()\`).

## Value

Integer vector of positions with an empty neighbour set.
