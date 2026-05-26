#' Effective sample size under spatial autocorrelation
#'
#' @description
#' Heuristic function to adjust sample size `n` by a spatial autocorrelation parameter `rho`.
#' Based on: Griffith, D.A., & Peres-Neto, P.R. (2006). Spatial Modeling in Ecology: The Flexibility of
#' Eigenfunction Spatial Analyses. Ecology, 87(10), 2603–2613.
#'
#' @param n integer, nominal sample size.
#' @param rho numeric in [0, 1], spatial autocorrelation intensity.
#' @return numeric, effective sample size.
#' @keywords internal
n_eff <- function(n, rho) {
  a <- 1 / (1 - exp(-1.92369))
  b <- (n - 1) / n
  c <- (1 - exp(-2.12373 * rho + 0.20024 * sqrt(rho)))
  n * (1 - a * b * c)
}

#' Get permutation runs
#' @param x ofemt_result
#' @return data.frame with columns Comparison, p_value, run, etc.
#' @keywords internal
get_perm_results <- function(x) {
  stopifnot(inherits(x, "ofemt_result"))
  x$perm_runs
}
