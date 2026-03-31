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
#' \dontrun{
#'   res <- ofemt(my_data, y = "Yield", x = "Treatment")
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
    cat(sprintf(
      "Cellsize: %s | Selected cells: %d\n",
      gi$Cellsize[1],
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
    cat("\n--- Means comparison ---\n")
    print(utils::head(x[["Means comparison"]], 10), row.names = FALSE)
    if (nrow(x[["Means comparison"]]) > 10) cat("... (truncated)\n")
  }

  if (!is.null(x[["ANOVA permutation test"]])) {
    cat(
      "\n--- Pairwise tests (median p-value across runs, corrected p-value) ---\n"
    )
    print(x[["ANOVA permutation test"]], row.names = FALSE)
  }

  invisible(x)
}
