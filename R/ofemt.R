#' OFE permutation analysis
#'
#' Analyzes unreplicated on-farm experiments to support field-specific
#' inference about treatment effects. Spatial statistical methods are
#' combined with permutation tests to compare the means of two or more
#' treatments.
#'
#' @param data An `sf` object containing point geometries and the response
#'   and treatment columns specified in `y` and `x`.
#' @param y response column (numeric).
#' @param x treatment column (factor/character).
#' @param cellsize,shift,angle_deg,buffer,min_per_cell grid settings (used when `grid` is NULL).
#' @param n_p number of permutations per ANOVA run.
#' @param n_s number of sampling runs.
#' @param alpha significance threshold for letters.
#' @param p_adjust_method p-value adjustment method across pairwise comparisons.
#'   The adjustment is applied *within each sampling run*, across the
#'   `choose(k, 2)` comparisons of that run; the reported `p_adj` is the median
#'   of those per-run adjusted values (see Details).
#' @param keep_components What spatial components to embed in the result. One of:
#'   \describe{
#'     \item{`"none"` (default)}{Nothing spatial is stored. Smallest object, but
#'       [plot_grid_selection()] and [plot.ofemt_result()] cannot be used on it.}
#'     \item{`"light"`}{Adds `grid`: the complete `ofe_grid` object used for the
#'       analysis (`grid_all`, `grid_sel`, `cell_stats`, `params`). Enough to
#'       redraw the full grid and the selected cells, but the observations
#'       themselves are *not* stored, so no points can be overlaid.}
#'     \item{`"full"`}{Everything in `"light"`, plus `cell_medians` (one point
#'       per selected cell carrying the per-cell median response, its treatment
#'       and the ANOVA residual used for the spatial diagnostics) and
#'       `points_joined` (every observation that fell inside a selected cell,
#'       with its `CellID`). This is the option that lets `plot()` overlay the
#'       raw points on the grid, which is the useful view when tuning
#'       `shift`, `angle_deg`, `buffer` or `cellsize`.}
#'   }
#'   Roughly, `"light"` costs one polygon per grid cell and `"full"` adds one
#'   row per observation.
#' @param grid optional: an `ofe_grid` (e.g., from [make_ofe_grid()]). If `NULL`,
#'   the grid is generated internally using `cellsize`, `min_per_cell`,
#'   `shift`, `angle_deg` and `buffer`, and an informative message is emitted.
#' @param crs optional target projected CRS if `data` is in lon/lat.
#' @param seed integer; controls reproducibility of the permutation sampling
#'   (the grid itself is deterministic from its construction arguments). Set
#'   `seed = NULL` to let results vary across runs.
#'
#' @export
#'
#' @details
#' The OFE-mean test accounts for spatial dependence when comparing
#' treatments in unreplicated on-farm experiments. The procedure involves:
#'
#' \enumerate{
#'   \item Aggregating the georeferenced observations within grid cells.
#'   \item Estimating the spatial autocorrelation of the treatment-adjusted
#'         residuals.
#'   \item Calculating the effective sample size (ESS) from the estimated
#'         spatial dependence.
#'   \item Drawing repeated balanced subsamples whose size is determined by
#'         the ESS.
#'   \item Performing pairwise permutation analysis of variance tests for
#'         each subsample.
#'   \item Generating an empirical distribution of p-values for each
#'         treatment comparison.
#' }
#'
#' The median of each empirical p-value distribution is reported as the
#' p-value associated with the null hypothesis of no treatment effect.
#'
#' For experiments with more than two treatments, all pairwise treatment
#' comparisons are performed.
#'
#' \strong{Multiplicity adjustment}
#'
#' When `p_adjust_method != "none"`, multiplicity adjustment is performed
#' separately within each sampling run. [stats::p.adjust()] is applied to the
#' `choose(k, 2)` pairwise p-values obtained in that run, where `k` is the
#' number of treatments. The reported adjusted p-value is the median of the
#' resulting empirical distribution of adjusted p-values.
#'
#' The `p_adj` column in `perm_runs` contains the adjusted p-value from each
#' sampling run. Therefore, [plot_pvalue_hist()] displays the empirical
#' distribution used to calculate the reported median.
#'
#' \strong{Compact letter display}
#'
#' Treatments are ordered by decreasing response before the compact letter
#' display is generated. Treatments that do not share a letter are considered
#' significantly different at the significance level specified by `alpha`.
#' Treatment labels are internally recoded before calling
#' [multcompView::multcompLetters()] and subsequently restored. Consequently,
#' labels containing spaces or special characters are preserved in the
#' results.
#'
#' @return An object of class `ofemt_result`: a list with
#'   \describe{
#'     \item{`General information`}{One-row data frame with the cell size, the
#'       number of cells in the full grid and in the selection, the
#'       min/median/max number of observations per selected cell, the number of
#'       cells entering the analysis (`n`), the effective sample size (`ESS`),
#'       the spatial autocorrelation estimate (`Rho`) and Moran's *I*.}
#'     \item{`Cells per treatment`}{`table` of selected cells per treatment.}
#'     \item{`ANOVA permutation test`}{One row per pairwise comparison, with the
#'       median p-value across runs (`p_value`) and its multiplicity-adjusted
#'       counterpart (`p_adj`).}
#'     \item{`Means comparison`}{One row per treatment, sorted by decreasing
#'       median response, with the compact letter display.}
#'     \item{`perm_runs`}{The raw per-run output: `n_s * choose(k, 2)` rows with
#'       columns `Trt_1`, `Trt_2`, `Comparison`, `p_value` (the permutation
#'       p-value of that comparison in that run), `p_adj` (the same value after
#'       adjusting within the run) and `run` (run index, `1:n_s`). This is the
#'       empirical p-value distribution the method is built on; it is what
#'       [plot_pvalue_hist()] draws and what you would use to inspect the
#'       run-to-run variability behind the reported medians.}
#'     \item{`params`}{The settings actually used (grid geometry, `n_p`, `n_s`,
#'       `alpha`, `p_adjust_method`, `seed`, and `grid_source`, which records
#'       whether the grid was supplied or built internally).}
#'     \item{`grid`, `cell_medians`, `points_joined`}{Optional spatial
#'       components; see `keep_components`.}
#'   }
#'
#' @references
#' Córdoba, M., Paccioretti, P. and Balzarini, M. (2025).
#' A new method to compare treatments in unreplicated on-farm experimentation.
#' \emph{Precision Agriculture}, 26, Article 4.
#' \doi{10.1007/s11119-024-10206-0}
#'
#' @seealso [make_ofe_grid()], [plot_grid_selection()], [plot_pvalue_hist()]
#'
#' @examples
#' \donttest{
#'   res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment",
#'                cellsize = 10, min_per_cell = 4, alpha = 0.05)
#'
#'   # Keep the geometries to inspect how the grid lands on the points
#'   res <- ofemt(ofe_f2, y = "Yield_tn", x = "Treatment",
#'                cellsize = 10, min_per_cell = 4,
#'                keep_components = "full")
#'   plot(res)
#' }
ofemt <- function(
  data,
  y,
  x,
  cellsize = 10,
  min_per_cell = 4,
  n_p = 1000,
  n_s = 200,
  alpha = 0.05,
  shift = c(0, 0),
  p_adjust_method = c("none", "bonferroni", "holm", "BH"),
  crs = NULL,
  keep_components = c("none", "light", "full"),
  grid = NULL,
  angle_deg = 0,
  buffer = 0,
  seed = 7L
) {
  keep_components <- match.arg(
    if (is.logical(keep_components)) {
      if (keep_components) "light" else "none"
    } else {
      keep_components
    },
    c("none", "light", "full")
  )
  p_adjust_method <- match.arg(
    if (is.logical(p_adjust_method)) {
      if (p_adjust_method) "bonferroni" else "none"
    } else {
      p_adjust_method
    },
    c("none", "bonferroni", "holm", "BH")
  )

  stopifnot(inherits(data, "sf"))
  stopifnot(
    length(y) == 1,
    length(x) == 1,
    y %in% names(data),
    x %in% names(data)
  )
  if (!is.numeric(data[[y]])) {
    stop("`y` must be numeric.")
  }
  if (!(is.character(data[[x]]) || is.factor(data[[x]]))) {
    stop("`x` must be factor/character.")
  }

  if (sf::st_is_longlat(data)) {
    if (is.null(crs)) {
      stop("`data` is lon/lat. Provide a projected `crs`.")
    }
    crs_obj <- sf::st_crs(crs)
    if (sf::st_is_longlat(crs_obj)) {
      stop("Target `crs` must be projected.")
    }
    data <- sf::st_transform(data, crs_obj)
  }

  # --- Treatment labels --------------------------------------------------
  # Everything downstream (formulas, permuco, multcompView) works on opaque
  # placeholder ids so that labels containing spaces, "+", "-", accents or any
  # other special character survive untouched. Labels are restored on the way
  # out via `trt_label()`.
  trt_raw <- as.character(data[[x]])
  if (anyNA(trt_raw)) {
    stop("`x` contains missing values.")
  }
  trt_labels <- if (is.factor(data[[x]])) {
    lv <- levels(data[[x]])
    lv[lv %in% trt_raw]
  } else {
    sort(unique(trt_raw))
  }
  if (length(trt_labels) < 2) {
    stop("Need at least two treatments.")
  }
  trt_ids <- paste0("t", seq_along(trt_labels))
  trt_label <- function(id) trt_labels[match(id, trt_ids)]
  data[[".trt"]] <- trt_ids[match(trt_raw, trt_labels)]

  grid_source <- "generated"
  used_params <- list(
    cellsize = if (length(cellsize) == 1) {
      c(cellsize, cellsize)
    } else {
      cellsize[1:2]
    },
    shift = shift,
    angle_deg = angle_deg,
    buffer = buffer,
    min_per_cell = min_per_cell
  )
  grid_obj <- NULL

  # --- grid handling

  if (!is.null(grid) && !inherits(grid, "ofe_grid")) {
    stop("`grid` must be ofe_grid.")
  }

  if (!is.null(grid)) {
    stopifnot(inherits(grid$grid_sel, "sf"))
    grid_source <- "provided_ofe_grid_sel"
    grid_obj <- grid
    used_params <- grid$params
    grid <- grid$grid_sel
  }

  # --- generated path: call make_ofe_grid and DO NOT re-apply min threshold here
  if (grid_source == "generated") {
    message(
      "`grid` not provided: building one internally via `make_ofe_grid()`. ",
      "Pass a pre-built `ofe_grid` to inspect or reuse the selection."
    )
    grid_obj <- make_ofe_grid(
      data = data,
      x = x,
      cellsize = used_params$cellsize,
      min_per_cell = used_params$min_per_cell,
      angle_deg = used_params$angle_deg,
      buffer = used_params$buffer,
      shift = used_params$shift
    )
    grid <- grid_obj$grid_sel
    used_params <- grid_obj$params
  } else {
    if (is.na(sf::st_crs(grid))) {
      stop("CRS mismatch between `grid` and `data`.")
    }
    if (!any(sf::st_geometry_type(grid) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("`grid` must be polygons.")
    }
    if (sf::st_crs(grid) != sf::st_crs(data)) {
      warning("`grid` CRS was transformed to `data` CRS.", call. = FALSE)
      grid <- sf::st_transform(grid, sf::st_crs(data))
      grid_obj$grid_sel <- grid
      if (!is.null(grid_obj$grid_all)) {
        grid_obj$grid_all <- sf::st_transform(
          grid_obj$grid_all,
          sf::st_crs(data)
        )
      }
    }
  }
  min_used <- used_params$min_per_cell

  if (!"CellID" %in% names(grid)) {
    grid$CellID <- seq_len(nrow(grid))
  }

  # Join points to selected cells. `grid` here already contains only
  # single-treatment cells with >= min_per_cell observations (filtered by
  # select_grid()), so we trust grid_sel and do not re-filter.
  jdat <- sf::st_join(
    data,
    grid[, "CellID"],
    join = sf::st_intersects,
    left = FALSE
  )

  if (nrow(jdat) == 0L) {
    stop("No data points fall in the selected grid cells.")
  }

  # Cell medians (per-cell median yield + coordinates, keep treatment label).
  # Coordinates go into dot-prefixed columns so that user columns literally
  # named "X"/"Y" are not clobbered, and st_drop_geometry() is used directly to
  # avoid data.frame()'s check.names mangling of the response column name.
  coords <- sf::st_coordinates(jdat)
  jdat_df <- sf::st_drop_geometry(jdat)
  jdat_df[[".X"]] <- coords[, "X"]
  jdat_df[[".Y"]] <- coords[, "Y"]

  cell_keys <- list(CellID = jdat_df$CellID, .trt = jdat_df$.trt)
  cell_med <- stats::aggregate(
    jdat_df[, c(y, ".X", ".Y")],
    by = cell_keys,
    FUN = stats::median
  )
  cell_n <- stats::aggregate(
    jdat_df[[y]],
    by = cell_keys,
    FUN = length
  )
  names(cell_n)[3] <- "n_obs"
  cell_med <- merge(cell_med, cell_n, by = c("CellID", ".trt"))

  # Restore the original treatment column label (the one passed in `x`)
  cell_med[[x]] <- trt_label(cell_med$.trt)

  # One-way ANOVA on cell medians to get residuals for spatial diagnostics
  my_model <- stats::lm(
    stats::as.formula(paste(quote_name(y), "~", ".trt")),
    data = cell_med
  )
  cell_med$residuos <- stats::residuals(my_model)

  cell_sf <- sf::st_as_sf(
    cell_med,
    coords = c(".X", ".Y"),
    crs = sf::st_crs(data)
  )

  # Spatial weights: nearest neighbours, distance-weighted.
  # `spdep` warns whenever the resulting graph is not fully connected
  # ("neighbour object has N sub-graphs"). That is expected for OFE trial
  # shapes and carries no consequence for the diagnostics below, so it is
  # muffled; the condition that *does* matter — cells left without any
  # neighbour — is checked explicitly right after.
  k1 <- without_subgraph_warnings(
    spdep::knn2nb(spdep::knearneigh(cell_sf))
  )
  dmax <- max(unlist(spdep::nbdists(k1, cell_sf)))
  # `dnearneigh()`'s upper bound is exclusive (d1 <= d < d2). `dmax` is by
  # construction the nearest-neighbour distance of *some* cell, so passing it
  # verbatim drops exactly that cell from its own graph and leaves it with an
  # empty neighbour set. The relative nudge is ~7 orders of magnitude above
  # double precision and geometrically negligible.
  gri <- without_subgraph_warnings(
    spdep::dnearneigh(cell_sf, 0, dmax * (1 + 1e-9))
  )

  isolated <- nb_no_links(gri)
  if (length(isolated) == length(gri)) {
    stop(
      "No cell has a neighbour within ",
      format(dmax, digits = 4),
      " map units, so the spatial weights matrix is empty. ",
      "Increase `cellsize` or lower `min_per_cell` so that the selected ",
      "cells are closer together.",
      call. = FALSE
    )
  }
  zero_policy <- length(isolated) > 0L
  if (zero_policy) {
    warning(
      length(isolated),
      " of ",
      length(gri),
      " selected cells have no neighbour within ",
      format(dmax, digits = 4),
      " map units; they contribute nothing to Rho and Moran's I.",
      call. = FALSE
    )
  }

  dist <- spdep::nbdists(gri, cell_sf)
  fdist <- lapply(dist, function(d) 1 / (d / dmax))
  lw <- spdep::nb2listw(
    gri,
    glist = fdist,
    style = "W",
    zero.policy = zero_policy
  )

  # Spatial autocorrelation of residuals + effective sample size
  rho <- spatialreg::aple(cell_sf$residuos, lw)
  rho <- min(1, max(0, rho))
  moran_test <- spdep::moran.test(
    cell_sf$residuos,
    lw,
    zero.policy = zero_policy
  )
  MI <- unname(moran_test$estimate["Moran I statistic"])
  n <- nrow(cell_sf)
  ess <- round(n_eff(n = n, rho = rho), 0)

  # Treatment medians (on cell medians), ordered by decreasing response so that
  # the compact letter display runs "a" (highest) downwards.
  trt_med_df <- stats::aggregate(
    cell_med[[y]],
    by = list(.trt = cell_med$.trt),
    FUN = stats::median
  )
  names(trt_med_df) <- c(".trt", y)
  trt_med_df <- trt_med_df[
    order(trt_med_df[[y]], decreasing = TRUE),
    ,
    drop = FALSE
  ]

  trt_order <- trt_med_df$.trt
  n_trt <- length(trt_order)
  compar <- utils::combn(trt_order, 2, simplify = TRUE)
  comp_key <- data.frame(
    .trt1 = compar[1, ],
    .trt2 = compar[2, ],
    Trt_1 = trt_label(compar[1, ]),
    Trt_2 = trt_label(compar[2, ]),
    stringsAsFactors = FALSE
  )
  comp_key$Comparison <- paste(comp_key$Trt_1, "vs.", comp_key$Trt_2)
  comp_key$.id <- paste(comp_key$.trt1, comp_key$.trt2, sep = "-")

  perm_formula <- stats::as.formula(paste(quote_name(y), "~ .trt"))

  # Multiple-permutation runs
  multipermutacion <- function(p) {
    # Sample ceiling(ess / n_trt) cells per treatment without replacement
    base_perm_idx <- unlist(lapply(trt_order, function(t) {
      sample_values(which(cell_sf$.trt == t), ceiling(ess / n_trt))
    }))
    base_perm <- cell_sf[base_perm_idx, ]

    permt_trat <- function(j) {
      pair <- compar[, j]
      sub <- base_perm[base_perm$.trt %in% pair, ]
      perm_res <- suppressWarnings(
        permuco::aovperm(perm_formula, data = sub, np = n_p)
      )
      data.frame(
        .id = paste(pair, collapse = "-"),
        p_value = perm_res$table$`resampled P(>F)`[1],
        stringsAsFactors = FALSE
      )
    }

    out <- do.call(rbind, lapply(seq_len(ncol(compar)), permt_trat))
    # Adjust for multiplicity *within* the run, so the reported p_adj is the
    # median of a genuine distribution of adjusted p-values.
    out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust_method)
    out$run <- p
    out
  }

  perm_runs <- if (is.null(seed)) {
    do.call(rbind, lapply(seq_len(n_s), multipermutacion))
  } else {
    withr::with_seed(
      seed,
      do.call(rbind, lapply(seq_len(n_s), multipermutacion))
    )
  }

  # Median p-value per comparison (raw and adjusted), in comp_key order
  med_by <- function(v) {
    unname(tapply(
      v,
      factor(perm_runs$.id, levels = comp_key$.id),
      function(z) stats::median(z, na.rm = TRUE)
    ))
  }
  pvals_by_comp <- data.frame(
    Comparison = comp_key$Comparison,
    p_value = med_by(perm_runs$p_value),
    p_adj = med_by(perm_runs$p_adj),
    stringsAsFactors = FALSE
  )

  # User-facing perm_runs: real labels, no internal ids
  m <- match(perm_runs$.id, comp_key$.id)
  perm_runs <- data.frame(
    Trt_1 = comp_key$Trt_1[m],
    Trt_2 = comp_key$Trt_2[m],
    Comparison = comp_key$Comparison[m],
    p_value = perm_runs$p_value,
    p_adj = perm_runs$p_adj,
    run = perm_runs$run,
    stringsAsFactors = FALSE
  )

  # Compact letter display from adjusted p-values. Placeholder ids keep
  # multcompView's "-"-separated name parsing safe for any treatment label,
  # and `trt_order` (decreasing means) fixes the letter assignment order.
  pvec <- pvals_by_comp$p_adj
  names(pvec) <- comp_key$.id
  letras_comp <- multcompView::multcompLetters(pvec, threshold = alpha)
  my_letters <- if (!is.null(letras_comp$monospacedLetters)) {
    letras_comp$monospacedLetters
  } else {
    letras_comp$Letters
  }

  means_table <- trt_med_df
  means_table$letters <- unname(my_letters[match(
    means_table$.trt,
    names(my_letters)
  )])
  means_table$.trt <- trt_label(means_table$.trt)
  # Rename to the user-facing column names
  names(means_table)[names(means_table) == ".trt"] <- x
  names(means_table)[names(means_table) == y] <- paste0(y, "_mean")
  row.names(means_table) <- NULL

  cells_per_trt <- table(factor(cell_med[[x]], levels = trt_labels))

  # Use the selected-cell n_obs as the per-cell observation counts
  n_obs_sel <- cell_med$n_obs
  infogral <- data.frame(
    Cellsize = paste(used_params$cellsize, collapse = " x "),
    Total.Cells = if (is.null(grid_obj$grid_all)) {
      NA_integer_
    } else {
      nrow(grid_obj$grid_all)
    },
    Selected.Cells = nrow(cell_sf),
    Min.Obs.Cell = min(n_obs_sel),
    Median.Obs.Cell = stats::median(n_obs_sel),
    Max.Obs.Cell = max(n_obs_sel),
    n = n,
    ESS = ess,
    Rho = rho,
    MoranI = MI,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  result <- list(
    `General information` = infogral,
    `Cells per treatment` = cells_per_trt,
    `ANOVA permutation test` = pvals_by_comp,
    `Means comparison` = means_table,
    perm_runs = perm_runs,
    params = list(
      cellsize = used_params$cellsize,
      shift = used_params$shift,
      angle_deg = used_params$angle_deg,
      buffer = used_params$buffer,
      min_per_cell = min_used,
      n_p = n_p,
      n_s = n_s,
      alpha = alpha,
      p_adjust_method = p_adjust_method,
      seed = seed,
      grid_source = grid_source,
      keep_components = keep_components
    )
  )

  if (keep_components != "none") {
    result$grid <- grid_obj
    if (keep_components == "full") {
      result$cell_medians <- cell_sf
      result$points_joined <- jdat
    }
  }

  class(result) <- c("ofemt_result", "list")
  result
}
