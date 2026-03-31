#' Plot grid selection (base R)
#' @param grid_all sf with all grid cells (polygons).
#' @param grid_sel sf with selected cells. Optional.
#' @param pts sf points. Optional.
#' @export
plot_grid_selection <- function(
  ofe_grid #,
  # grid_all,
  # grid_sel = NULL,
  # pts = NULL
) {
  stopifnot(inherits(grid_all, "ofe_grid"))
  grid_all <- grid[['grid_all']]
  grid_sel <- grid[['grid_sel']]
  stopifnot(inherits(grid_all, "sf"))
  plot(sf::st_geometry(grid_all), border = "grey85", lwd = 1)
  if (nrow(grid_sel) > 0) {
    plot(sf::st_geometry(grid_sel), add = TRUE, border = "black", lwd = 1.5)
  }
  if (!is.null(pts)) {
    plot(sf::st_geometry(pts), add = TRUE, pch = 16, cex = .6)
  }
  legend(
    "topleft",
    legend = c("Full grid", "Selected cells", "Points"),
    lty = c(1, 1, NA),
    pch = c(NA, NA, 16),
    col = c("grey60", "black", "black"),
    bty = "n"
  )
}

#' Plot histogram(s) of permutation p-values per comparison
#'
#' This function produces one or several histograms showing the distribution of
#' permutation-derived *p*-values for each comparison included in an
#' [`ofemt_result`] object. A vertical dashed red line indicates the
#' significance threshold (`alpha`) used during OFE permutation analysis.
#'
#' @param results An object of class [`ofemt_result`], typically obtained from
#'   [ofemt()]. It must contain the elements `perm_runs` (a data frame with
#'   columns `"Comparison"` and `"p_value"`) and `params$alpha`.
#'
#' @return A [`ggplot2`] object representing one or several histograms of
#'   permutation *p*-values, faceted by comparison, with the alpha threshold
#'   indicated by a dashed red vertical line.
#'
#' @details
#' Each facet corresponds to a different comparison contained in
#' `results$perm_runs`. The function automatically adds a legend entry for the
#' alpha threshold (e.g., `"U+03b1 = 0.05"`) below the comparison legend.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'res' is an ofemt_result obtained from ofemt():
#'   plot_pvalue_hist(res, alpha = 0.05)
#' }
#'
#' @export
plot_pvalue_hist <- function(results) {
  if (!inherits(x, "ofemt_result")) {
    stop('results must be an ofemt_result object')
  }

  permitation_runs <- results$perm_runs
  params_alpha <- results$params$alpha

  stopifnot(
    is.data.frame(permitation_runs),
    all(c("Comparison", "p_value") %in% names(permitation_runs))
  )

  alpha_df <- data.frame(
    x = params_alpha,
    label = paste0("\\u03b1 = ", params_alpha)
  )

  ggplot2::ggplot(
    permitation_runs,
    ggplot2::aes(x = p_value, fill = Comparison, colour = Comparison)
  ) +
    ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
    ggplot2::facet_wrap(~Comparison, scales = "free_y") +
    ggplot2::geom_vline(
      data = alpha_df,
      ggplot2::aes(xintercept = x, linetype = label),
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
