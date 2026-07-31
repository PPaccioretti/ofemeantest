# Quote a column name for use in a formula

Column names coming from real field data routinely contain spaces,
\`+\`, accents or parentheses. Backticking them keeps
\[stats::as.formula()\] happy without forcing the user to rename
anything.

## Usage

``` r
quote_name(nm)
```

## Arguments

- nm:

  Character scalar, a column name.

## Value

Character scalar, the backtick-quoted name.
