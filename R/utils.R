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

#' Evaluate an expression muffling `spdep`'s connectivity notes
#'
#' @description
#' `spdep` emits a warning whenever a neighbour object is not fully connected
#' (`"neighbour object has N sub-graphs"`). For the cell-median point pattern
#' used here that is expected and harmless: the weights matrix is built from
#' every cell's nearest-neighbour distance, so disconnected components simply
#' reflect the shape of the trial. Only the *absence* of neighbours is a real
#' problem, and that is checked explicitly by the caller (see [nb_no_links()]).
#'
#' Warnings that do not mention sub-graphs are left untouched.
#'
#' @param expr Expression to evaluate.
#' @return The value of `expr`.
#' @keywords internal
without_subgraph_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("sub-?graph", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Regions with no neighbours in an `nb` object
#'
#' In `spdep`, a region with an empty neighbour set is stored as the integer
#' `0L`. This helper returns the positions of those regions.
#'
#' @param nb An `nb` object (e.g. from `spdep::dnearneigh()`).
#' @return Integer vector of positions with an empty neighbour set.
#' @keywords internal
nb_no_links <- function(nb) {
  which(vapply(
    nb,
    function(z) length(z) == 1L && !is.na(z[1L]) && z[1L] == 0L,
    logical(1)
  ))
}

#' Quote a column name for use in a formula
#'
#' Column names coming from real field data routinely contain spaces, `+`,
#' accents or parentheses. Backticking them keeps [stats::as.formula()] happy
#' without forcing the user to rename anything.
#'
#' @param nm Character scalar, a column name.
#' @return Character scalar, the backtick-quoted name.
#' @keywords internal
quote_name <- function(nm) {
  paste0("`", gsub("`", "\\\\`", nm), "`")
}

#' Sample without the `sample()` length-one trap
#'
#' `sample(x, k)` treats a length-one numeric `x` as `1:x`. This helper always
#' samples *from the elements* of `x`.
#'
#' @param x Vector to sample from.
#' @param size Number of elements to draw.
#' @return A vector of `size` elements of `x`.
#' @keywords internal
sample_values <- function(x, size) {
  x[sample.int(length(x), size = min(size, length(x)))]
}
