# On-Farm single strip treatment trial (companion dataset)

Companion on-farm experiment used to demonstrate the package on a second
field. Same layout as \[ofe_f2\]: a single fertilized strip compared
against an adjacent control strip on dense yield-monitor data, with
yields cleaned following Vega et al. (2019).

## Usage

``` r
ofe_p
```

## Format

A sf object with 3070 rows and 3 variables:

- Treatment:

  character, identifying Fertilized or Control strip

- Yield_tn:

  numeric, grain yield in tn/ha

- geom:

  sf geometry column

## Details

Coordinate reference system is "WGS 84 / UTM zone 20S", EPSG:32720.
