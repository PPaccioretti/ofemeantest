#' On-Farm single strip treatment trial (Field F2)
#'
#' The impact of broadcast phosphorus (P) fertilization on corn was evaluated.
#' The Fertilized strip received superphosphate at rates ranging from 250 to
#' 500 kg ha-1 \[0-21-0\], while the Control strip received no P fertilizer
#' (0 kg ha-1). Each strip was 2.2 ha in area, totaling 100 ha.
#'
#' @format A sf object with 3070 rows and 3 variables:
#' \describe{
#'   \item{Treatment}{character, identifying Fertilized or Control strip}
#'   \item{Yield_tn}{numeric, corn grain yield in tn/ha}
#'   \item{geom}{sf geometry column}
#' }
#' @details
#' Coordinate reference system is "WGS 84 / UTM zone 20S", EPSG:32720.
"ofe_f2"

#' On-Farm single strip treatment trial (companion dataset)
#'
#' Companion on-farm experiment used to demonstrate the package on a second
#' field. Same layout as [ofe_f2]: a single fertilized strip compared against
#' an adjacent control strip on dense yield-monitor data, with yields cleaned
#' following Vega et al. (2019).
#'
#' @format A sf object with 3070 rows and 3 variables:
#' \describe{
#'   \item{Treatment}{character, identifying Fertilized or Control strip}
#'   \item{Yield_tn}{numeric, grain yield in tn/ha}
#'   \item{geom}{sf geometry column}
#' }
#' @details
#' Coordinate reference system is "WGS 84 / UTM zone 20S", EPSG:32720.
"ofe_p"
