#' Plot grid selection (base R)
#'
#' Quick visual check of an `ofe_grid` object: full grid in grey, selected cells
#' (those that passed the single-treatment and `min_per_cell` filters) outlined
#' in black. Optionally overlays the point data.
#'
#' @param ofe_grid An object of class `ofe_grid` (typically returned by
#'   [make_ofe_grid()]).
#' @param data Optional `sf` points overlaid on the grid.
#'
#' @return Invisibly returns `NULL`. Called for its side effect (a base R plot).
#' @export
plot_grid_selection <- function(ofe_grid, data = NULL) {
  stopifnot(inherits(ofe_grid, "ofe_grid"))

  grid_all <- ofe_grid[["grid_all"]]
  grid_sel <- ofe_grid[["grid_sel"]]
  stopifnot(inherits(grid_all, "sf"))

  plot(sf::st_geometry(grid_all), border = "grey85", lwd = 1)
  if (!is.null(grid_sel) && nrow(grid_sel) > 0) {
    plot(sf::st_geometry(grid_sel), add = TRUE, border = "black", lwd = 1.5)
  }
  if (!is.null(data)) {
    stopifnot(inherits(data, "sf"))
    plot(sf::st_geometry(data), add = TRUE, pch = 16, cex = 0.6)
  }
  legend(
    "topleft",
    legend = c("Full grid", "Selected cells", "Points"),
    lty = c(1, 1, NA),
    pch = c(NA, NA, 16),
    col = c("grey60", "black", "black"),
    bty = "n"
  )
  invisible(NULL)
}

#' Plot histogram(s) of permutation p-values per comparison
#'
#' Produces one or several histograms showing the distribution of
#' permutation-derived *p*-values for each comparison contained in an
#' [`ofemt_result`] object. A dashed red vertical line marks the significance
#' threshold (`alpha`) used in the analysis.
#'
#' @param results An object of class `ofemt_result`, typically obtained from
#'   [ofemt()]. Must contain `perm_runs` (data frame with `Comparison` and
#'   `p_value`) and `params$alpha`.
#'
#' @return A `ggplot` object if **ggplot2** is installed; otherwise the
#'   function falls back to base R `hist()` and returns `NULL` invisibly.
#'
#' @examples
#' \dontrun{
#'   res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment", cellsize = 9)
#'   plot_pvalue_hist(res)
#' }
#'
#' @export
plot_pvalue_hist <- function(results) {
  # Declared here to silence R CMD check's "no visible binding for global
  # variable '.data'" note without adding rlang as an Imports dependency.
  # ggplot2 still resolves `.data$col` correctly at evaluation time via its
  # data mask, which shadows this local binding.
  .data <- NULL

  if (!inherits(results, "ofemt_result")) {
    stop("`results` must be an ofemt_result object")
  }

  perm_runs <- results$perm_runs
  alpha <- results$params$alpha

  stopifnot(
    is.data.frame(perm_runs),
    all(c("Comparison", "p_value") %in% names(perm_runs))
  )

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback: one panel per comparison
    comps <- unique(perm_runs$Comparison)
    op <- graphics::par(mfrow = c(1, length(comps)))
    on.exit(graphics::par(op), add = TRUE)
    for (cmp in comps) {
      graphics::hist(
        perm_runs$p_value[perm_runs$Comparison == cmp],
        main = cmp,
        xlab = "p-value",
        col = "grey80",
        border = "white"
      )
      graphics::abline(v = alpha, lty = 2, col = "red")
    }
    return(invisible(NULL))
  }

  alpha_df <- data.frame(
    x = alpha,
    label = paste0("α = ", alpha)
  )

  ggplot2::ggplot(
    perm_runs,
    ggplot2::aes(x = .data$p_value, fill = .data$Comparison, colour = .data$Comparison)
  ) +
    ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
    ggplot2::facet_wrap(~ .data$Comparison, scales = "free_y") +
    ggplot2::geom_vline(
      data = alpha_df,
      ggplot2::aes(xintercept = .data$x, linetype = .data$label),
      color = "red"
    ) +
    ggplot2::scale_linetype_manual(values = "dashed", name = NULL) +
    ggplot2::labs(y = "Absolute frequency", x = "p-value") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 2)
    )
}
