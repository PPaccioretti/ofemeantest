# Sample without the \`sample()\` length-one trap

\`sample(x, k)\` treats a length-one numeric \`x\` as \`1:x\`. This
helper always samples \*from the elements\* of \`x\`.

## Usage

``` r
sample_values(x, size)
```

## Arguments

- x:

  Vector to sample from.

- size:

  Number of elements to draw.

## Value

A vector of \`size\` elements of \`x\`.
