#' Plot the grid, the selected cells and the observations
#'
#' Draws the three layers that matter when tuning a grid, on one set of axes:
#' the **full grid** in light grey, the **selected cells** (those that passed
#' the single-treatment and `min_per_cell` filters) shaded and outlined in
#' black, and the **observations** as points. Seeing the points against the
#' cell boundaries is the quickest way to judge the effect of `cellsize`,
#' `shift`, `angle_deg` and `buffer`.
#'
#' @param x Either an `ofe_grid` (from [make_ofe_grid()]) or an `ofemt_result`
#'   (from [ofemt()]) that was run with `keep_components = "light"` or
#'   `"full"`. With an `ofemt_result` every layer is taken from the object
#'   itself, so no extra arguments are needed; the points are only available
#'   under `keep_components = "full"`.
#' @param data Optional `sf` points to overlay. Rarely needed: the observations
#'   are taken from the object itself — `points_sel` for an `ofe_grid` built
#'   with `make_ofe_grid(return_points = TRUE)`, `points_joined` for an
#'   `ofemt_result` run with `keep_components = "full"`. Pass `data` when the
#'   object carries no points, or to override the stored ones.
#' @param points Logical; set to `FALSE` to skip the point layer. When `TRUE`
#'   (the default) and no observations are available, a message explains how to
#'   obtain them rather than silently drawing a grid without points.
#' @param main Plot title. Defaults to a one-line summary of the grid
#'   parameters, which is what makes successive calls comparable.
#' @param legend Logical; draw the legend. Default `TRUE`.
#' @param legend_pos Where to place the legend, passed to [graphics::legend()]
#'   (e.g. `"topleft"`, `"bottomright"`, or `"top"`). Default `"topleft"`.
#' @param ... Further arguments passed to the underlying [plot()] call for the
#'   full-grid layer.
#'
#' @return Invisibly returns `NULL`. Called for its side effect (a base R plot).
#'
#' @examples
#' \dontrun{
#'   g <- make_ofe_grid(ofe_f2, x = "Treatment", cellsize = 9, min_per_cell = 4)
#'   plot_grid_selection(g, data = ofe_f2)
#'   plot(g, data = ofe_f2)          # same thing
#'
#'   res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment", cellsize = 9,
#'                keep_components = "full")
#'   plot(res)                        # grid + selection + points, no extra args
#' }
#'
#' @seealso [make_ofe_grid()], [ofemt()]
#' @importFrom graphics legend
#' @importFrom grDevices adjustcolor
#' @export
plot_grid_selection <- function(
  x,
  data = NULL,
  points = TRUE,
  main = NULL,
  legend = TRUE,
  legend_pos = "topleft",
  ...
) {
  if (inherits(x, "ofemt_result")) {
    grid_obj <- x[["grid"]]
    if (is.null(grid_obj)) {
      stop(
        "This `ofemt_result` carries no geometries. Re-run `ofemt()` with ",
        "`keep_components = \"full\"` (grid + points) or `\"light\"` (grid only).",
        call. = FALSE
      )
    }
    if (is.null(data) && isTRUE(points)) {
      data <- x[["points_joined"]]
    }
    no_points_hint <- paste0(
      "re-run `ofemt()` with `keep_components = \"full\"`, ",
      "or pass the observations via `data =`"
    )
    if (is.null(main)) {
      main <- grid_params_label(x[["params"]])
    }
  } else if (inherits(x, "ofe_grid")) {
    grid_obj <- x
    # `make_ofe_grid(return_points = TRUE)` stores the observations that fell
    # in the selected cells as `points_sel`; use them so the argument actually
    # has a visible effect here.
    if (is.null(data) && isTRUE(points)) {
      data <- x[["points_sel"]]
    }
    no_points_hint <- paste0(
      "rebuild the grid with `make_ofe_grid(..., return_points = TRUE)`, ",
      "or pass the observations via `data =`"
    )
    if (is.null(main)) {
      main <- grid_params_label(x[["params"]])
    }
  } else {
    stop("`x` must be an `ofe_grid` or an `ofemt_result`.", call. = FALSE)
  }

  grid_sel <- grid_obj[["grid_sel"]]
  grid_all <- grid_obj[["grid_all"]] %||% grid_sel
  stopifnot(inherits(grid_all, "sf"))

  if (!isTRUE(points)) {
    data <- NULL
  } else if (is.null(data)) {
    # Points were asked for (the default) but the object carries none. Say so
    # instead of silently drawing a grid with no observations on it.
    message(
      "`plot_grid_selection()`: no observations to draw; ",
      no_points_hint,
      "."
    )
  }

  plot(
    sf::st_geometry(grid_all),
    border = "grey85",
    lwd = 0.7,
    main = main,
    cex.main = 0.85,
    font.main = 1,
    ...
  )
  has_sel <- !is.null(grid_sel) && nrow(grid_sel) > 0
  if (has_sel) {
    plot(
      sf::st_geometry(grid_sel),
      add = TRUE,
      border = "black",
      col = adjustcolor("#1f78b4", alpha.f = 0.20),
      lwd = 1.2
    )
  }
  has_pts <- !is.null(data)
  if (has_pts) {
    stopifnot(inherits(data, "sf"))
    plot(
      sf::st_geometry(data),
      add = TRUE,
      pch = 16,
      cex = 0.35,
      col = adjustcolor("black", alpha.f = 0.55)
    )
  }

  if (isTRUE(legend)) {
    # Filled squares mirror what is actually drawn, so the mapping needs no
    # explanation. The box is opaque: with `bty = "n"` the labels landed on
    # top of the grid and became unreadable.
    keep <- c(TRUE, has_sel, has_pts)
    graphics::legend(
      legend_pos,
      legend = c("Unselected cells", "Selected cells", "Observations")[keep],
      pch = c(22, 22, 16)[keep],
      pt.bg = c(
        "white",
        adjustcolor("#1f78b4", alpha.f = 0.20),
        NA
      )[keep],
      pt.cex = c(1.6, 1.6, 0.9)[keep],
      col = c("grey85", "black", adjustcolor("black", alpha.f = 0.55))[keep],
      bty = "o",
      bg = "white",
      box.col = "grey70",
      box.lwd = 0.8,
      cex = 0.8,
      inset = 0.01
    )
  }
  invisible(NULL)
}

#' @rdname plot_grid_selection
#' @export
plot.ofe_grid <- function(x, ...) {
  plot_grid_selection(x, ...)
}

#' @rdname plot_grid_selection
#' @export
plot.ofemt_result <- function(x, ...) {
  plot_grid_selection(x, ...)
}

#' One-line description of the grid parameters
#'
#' @param params The `params` list of an `ofe_grid` or `ofemt_result`.
#' @return A character scalar, or `NULL` when `params` is empty.
#' @keywords internal
grid_params_label <- function(params) {
  if (is.null(params) || !length(params)) {
    return(NULL)
  }
  num <- function(v) paste(format(v, trim = TRUE), collapse = " x ")
  bits <- c(
    if (!is.null(params$cellsize)) paste0("cellsize ", num(params$cellsize)),
    if (!is.null(params$shift)) paste0("shift ", num(params$shift)),
    if (!is.null(params$angle_deg)) {
      paste0("angle ", num(params$angle_deg), "\u00b0")
    },
    if (!is.null(params$buffer)) paste0("buffer ", num(params$buffer)),
    if (!is.null(params$min_per_cell)) {
      paste0("min/cell ", num(params$min_per_cell))
    }
  )
  if (!length(bits)) NULL else paste(bits, collapse = " | ")
}

#' Plot histogram(s) of permutation p-values per comparison
#'
#' Produces one histogram per pairwise comparison showing the empirical
#' distribution of permutation *p*-values across the `n_s` sampling runs of an
#' [`ofemt_result`]. Two reference lines are drawn on each panel: a solid line
#' at the **median** of the distribution — the value [ofemt()] reports for that
#' comparison — and a dashed line at the significance threshold **alpha**.
#'
#' @param results An object of class `ofemt_result`, typically obtained from
#'   [ofemt()]. Must contain `perm_runs` and `params$alpha`.
#' @param which Which p-values to plot. `"auto"` (the default) uses the
#'   multiplicity-adjusted values whenever the analysis was run with
#'   `p_adjust_method != "none"`, and the raw values otherwise. `"adjusted"`
#'   and `"raw"` force one or the other. The adjustment is applied within each
#'   run before the histogram is built, so the median line coincides with the
#'   `p_adj` reported in the `ANOVA permutation test` table.
#' @param bins Number of histogram bins. Default 30.
#'
#' @return A `ggplot` object if **ggplot2** is installed; otherwise the
#'   function falls back to base R `hist()` and returns `NULL` invisibly.
#'
#' @examples
#' \dontrun{
#'   res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment", cellsize = 9,
#'                p_adjust_method = "bonferroni")
#'   plot_pvalue_hist(res)              # adjusted p-values
#'   plot_pvalue_hist(res, which = "raw")
#' }
#'
#' @export
plot_pvalue_hist <- function(
  results,
  which = c("auto", "adjusted", "raw"),
  bins = 30
) {
  # Declared here to silence R CMD check's "no visible binding for global
  # variable '.data'" note without adding rlang as an Imports dependency.
  # ggplot2 still resolves `.data$col` correctly at evaluation time via its
  # data mask, which shadows this local binding.
  .data <- NULL

  if (!inherits(results, "ofemt_result")) {
    stop("`results` must be an ofemt_result object")
  }
  which <- match.arg(which)

  perm_runs <- results$perm_runs
  alpha <- results$params$alpha
  adjust_method <- results$params$p_adjust_method %||% "none"

  stopifnot(
    is.data.frame(perm_runs),
    all(c("Comparison", "p_value") %in% names(perm_runs))
  )

  use_adjusted <- switch(
    which,
    auto = !identical(adjust_method, "none") && "p_adj" %in% names(perm_runs),
    adjusted = TRUE,
    raw = FALSE
  )
  if (use_adjusted && !"p_adj" %in% names(perm_runs)) {
    stop(
      "`perm_runs` has no `p_adj` column; re-run `ofemt()` or use ",
      "`which = \"raw\"`.",
      call. = FALSE
    )
  }
  pcol <- if (use_adjusted) "p_adj" else "p_value"
  xlab <- if (use_adjusted) {
    paste0("Adjusted p-value (", adjust_method, ")")
  } else {
    "p-value"
  }

  perm_runs$.p <- perm_runs[[pcol]]

  medians <- stats::aggregate(
    perm_runs$.p,
    by = list(Comparison = perm_runs$Comparison),
    FUN = function(z) stats::median(z, na.rm = TRUE)
  )
  names(medians) <- c("Comparison", "x")

  median_lab <- "median p"
  alpha_lab <- paste0("\u03b1 = ", alpha)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback: one panel per comparison
    comps <- unique(perm_runs$Comparison)
    op <- graphics::par(mfrow = c(1, length(comps)))
    on.exit(graphics::par(op), add = TRUE)
    for (cmp in comps) {
      graphics::hist(
        perm_runs$.p[perm_runs$Comparison == cmp],
        main = cmp,
        xlab = xlab,
        col = "grey80",
        border = "white"
      )
      graphics::abline(
        v = medians$x[medians$Comparison == cmp],
        lty = 1,
        lwd = 2,
        col = "#1f78b4"
      )
      graphics::abline(v = alpha, lty = 2, lwd = 2, col = "red")
      graphics::legend(
        "topright",
        legend = c(median_lab, alpha_lab),
        lty = c(1, 2),
        col = c("#1f78b4", "red"),
        bty = "n",
        cex = 0.8
      )
    }
    return(invisible(NULL))
  }

  lines_df <- rbind(
    data.frame(
      Comparison = medians$Comparison,
      x = medians$x,
      label = median_lab,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Comparison = medians$Comparison,
      x = alpha,
      label = alpha_lab,
      stringsAsFactors = FALSE
    )
  )
  lines_df$label <- factor(lines_df$label, levels = c(median_lab, alpha_lab))

  ggplot2::ggplot(perm_runs, ggplot2::aes(x = .data$.p)) +
    ggplot2::geom_histogram(
      fill = "grey75",
      colour = "white",
      bins = bins
    ) +
    ggplot2::facet_wrap(~ .data$Comparison, scales = "free_y") +
    ggplot2::geom_vline(
      data = lines_df,
      ggplot2::aes(
        xintercept = .data$x,
        colour = .data$label,
        linetype = .data$label
      ),
      linewidth = 0.7
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(c("#1f78b4", "red"), c(median_lab, alpha_lab)),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(c("solid", "dashed"), c(median_lab, alpha_lab)),
      name = NULL
    ) +
    ggplot2::labs(y = "Absolute frequency", x = xlab)
}
