# Evaluate an expression muffling \`spdep\`'s connectivity notes

\`spdep\` emits a warning whenever a neighbour object is not fully
connected (\`"neighbour object has N sub-graphs"\`). For the cell-median
point pattern used here that is expected and harmless: the weights
matrix is built from every cell's nearest-neighbour distance, so
disconnected components simply reflect the shape of the trial. Only the
\*absence\* of neighbours is a real problem, and that is checked
explicitly by the caller (see \[nb_no_links()\]).

Warnings that do not mention sub-graphs are left untouched.

## Usage

``` r
without_subgraph_warnings(expr)
```

## Arguments

- expr:

  Expression to evaluate.

## Value

The value of \`expr\`.
