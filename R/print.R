#' Print method for `ofemt_result` objects
#'
#' Displays a concise summary of an [`ofemt_result`] object, including general
#' information about the experimental grid, descriptive statistics, and
#' permutation-based test results.
#'
#' @param x An object of class [`ofemt_result`], typically returned by
#'   [ofemt()]. The object must contain, at minimum, a `"General information"`
#'   data frame. Optionally, it may include components such as
#'   `"Means comparison"` and `"ANOVA permutation test"`.
#' @param ... Additional arguments passed to or from other methods (ignored).
#'
#' @details
#' This print method provides a human-readable summary in the console:
#' \itemize{
#'   \item **General information:** Cell size, number of selected cells, and
#'     per-cell observation counts (minimum, median, maximum).
#'   \item **Spatial statistics:** Sample size (`n`), effective sample size
#'     (`ESS`), spatial correlation (`Rho`), and Moran’s *I*.
#'   \item **Means comparison:** Shows the first 10 rows of the table of
#'     comparisons between treatment means.
#'   \item **ANOVA permutation test:** Displays pairwise tests summarizing
#'     median *p*-values and corrected *p*-values across permutation runs.
#' }
#'
#' If any of these components are missing, they are simply skipped in the
#' printed output.
#'
#' @return Invisibly returns the input object `x`, unchanged.
#'
#' @examples
#' \donttest{
#'   res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment")
#'   print(res)
#' }
#'
#' @export
print.ofemt_result <- function(x, ...) {
  if (!inherits(x, "ofemt_result")) {
    stop("x must be an ofemt_result object")
  }

  gi <- x[["General information"]]
  cat("\n=== OFE permutation analysis ===\n")

  if (is.data.frame(gi)) {
    total <- if ("Total.Cells" %in% names(gi)) gi$Total.Cells[1] else NA
    cat(sprintf(
      "Cellsize: %s | %sSelected cells: %d\n",
      gi$Cellsize[1],
      if (is.na(total)) "" else sprintf("Total cells: %d | ", total),
      gi$Selected.Cells[1]
    ))
    cat(sprintf(
      "Obs/cell (min/median/max): %s / %s / %s\n",
      gi$Min.Obs.Cell[1],
      gi$Median.Obs.Cell[1],
      gi$Max.Obs.Cell[1]
    ))
    cat(sprintf(
      "n: %d | ESS: %d | Rho: %.3f | Moran's I: %.3f\n",
      gi$n[1],
      gi$ESS[1],
      gi$Rho[1],
      gi$MoranI[1]
    ))
  }

  if (!is.null(x[["Means comparison"]])) {
    cat("\n--- Means comparison (sorted by decreasing mean) ---\n")
    print(utils::head(x[["Means comparison"]], 10), row.names = FALSE)
    if (nrow(x[["Means comparison"]]) > 10) cat("... (truncated)\n")
  }

  if (!is.null(x[["ANOVA permutation test"]])) {
    adj <- x[["params"]][["p_adjust_method"]]
    cat(sprintf(
      "\n--- Pairwise tests (median p across runs%s) ---\n",
      if (is.null(adj) || identical(adj, "none")) {
        ""
      } else {
        sprintf("; p_adj = median of per-run %s-adjusted p", adj)
      }
    ))
    print(x[["ANOVA permutation test"]], row.names = FALSE)
  }

  invisible(x)
}
