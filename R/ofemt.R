#' OFE permutation analysis
#'
#' Runs a cell-based OFE analysis with permutation ANOVA and spatial diagnostics.
#'
#' @param data sf points; must contain columns `y` and `x`.
#' @param y response column (numeric).
#' @param x treatment column (factor/character).
#' @param cellsize,shift,angle_deg,buffer,min_per_cell grid settings (used when `grid` is NULL).
#' @param n_p number of permutations per ANOVA run.
#' @param n_s number of sampling runs.
#' @param alpha significance threshold for letters.
#' @param p_adjust_method p-value adjustment method across pairwise comparisons.
#' @param keep_components `"none"`, `"light"`, or `"full"` to embed geometries in output.
#' @param grid optional: an `ofe_grid`, a list with `grid_sel`, or an `sf` object of polygons.
#' @param crs optional target projected CRS if `data` is in lon/lat.
#' @param nmin_cell,alpha_bonferroni **Deprecated.**
#'   Use `min_per_cell` instead of `nmin_cell`, and `p_adjust_method = "bonferroni"`
#'   instead of `alpha_bonferroni`.
#'
#' @return An object of class `ofemt_result`.
#' @export
#'
#' @details
#' Both `nmin_cell` and `alpha_bonferroni` are retained for backward compatibility
#' but will be removed in a future release.
#' Use `min_per_cell` to define the minimum number of observations per grid cell,
#' and control Bonferroni correction through the argument `p_adjust_method`.
#'
#' @examples
#' \dontrun{
#'   my_data <- ofe_f2
#'   res <- ofemt(my_data, y = "Yield", x = "Treatment",
#'                cellsize = 10, min_per_cell = 4, alpha = 0.05)
#' }
ofemt <- function(
  data,
  y,
  x,
  cellsize = 10,
  min_per_cell = 4,
  nmin_cell = NULL,
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
  alpha_bonferroni = NULL
) {
  # --- Deprecation handling ---
  if (!is.null(nmin_cell)) {
    warning(
      "`nmin_cell` is deprecated; use `min_per_cell` instead.",
      call. = FALSE
    )
    min_per_cell <- nmin_cell
  }

  if (!is.null(alpha_bonferroni)) {
    warning(
      "`alpha_bonferroni` is deprecated; use `p_adjust_method = 'bonferroni'` instead.",
      call. = FALSE
    )
    p_adjust_method <- "bonferroni"
  }

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

  data[[".trt"]] <- gsub(" ", ".", as.character(data[[x]]))
  treatments <- length(unique(data$.trt))
  if (treatments < 2) {
    stop("Need at least two treatments.")
  }

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
  # min_applied <- TRUE
  min_used <- min_per_cell
  info_note <- NA_character_

  # --- grid handling

  if (!is.null(grid) & !inherits(grid, "ofe_grid")) {
    stop("`grid` must be ofe_grid.")
  }

  if (!is.null(grid) & inherits(grid, "ofe_grid")) {
    stopifnot(inherits(grid$grid_sel, "sf"))
    grid_source <- "provided_ofe_grid_sel"
    used_params <- grid$params
    grid <- grid$grid_sel
    # class(grid) <- c('ofe_grid', class(grid))
    # min_applied <- FALSE
    min_used <- used_params$min_per_cell
  }

  # --- generated path: call make_ofe_grid and DO NOT re-apply min threshold here
  if (grid_source == "generated") {
    gint <- make_ofe_grid(
      data = data,
      x = x,
      cellsize = used_params$cellsize,
      angle_deg = used_params$angle_deg,
      buffer = used_params$buffer,
      shift = used_params$shift #,
      # min_per_cell = used_params$min_per_cell,
      # return_both = TRUE
    )
    grid <- gint$grid_sel
    ap <- attr(grid, "ofe_params", exact = TRUE)
    if (!is.null(ap)) {
      used_params <- ap
    }
    # min_applied <- FALSE
    min_used <- used_params$min_per_cell
  } else {
    if (is.na(sf::st_crs(grid))) {
      stop("CRS mismatch between `grid` and `data`.")
    }
    if (!any(sf::st_geometry_type(grid) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("`grid` must be polygons.")
    }
    if (sf::st_crs(grid) != sf::st_crs(data)) {
      warning("`grid` CRS was trasformed to `data` CRS.", call. = FALSE)
      grid <- sf::st_transform(grid, sf::st_crs(data))
    }
  }

  # if (!"CellID" %in% names(grid)) {
  #   grid$CellID <- seq_len(nrow(grid))
  # }

  # # join points to cells
  jdat <- sf::st_join(data, grid, join = sf::st_intersects, left = FALSE)

  # tmp <- stats::aggregate(
  #   jdat$.trt,
  #   by = list(CellID = jdat$CellID),
  #   FUN = function(z) {
  #     c(
  #       n_trt = length(unique(z)),
  #       trt = ifelse(length(unique(z)) == 1, unique(z), NA),
  #       n = length(z)
  #     )
  #   }
  # )
  # tmp <- do.call(data.frame, tmp)
  # names(tmp) <- c("CellID", "unique_trt", "trt", "n")

  # # eligibility: if min already applied upstream, only require 1 unique trt
  # eligible <- tmp$unique_trt == 1
  # # if (isFALSE(min_applied)) {
  # #   tmp$unique_trt == 1
  # # } else {
  # #   tmp$unique_trt == 1 & tmp$n >= min_per_cell
  # # }
  # if (!any(eligible)) {
  #   stop(
  #     "No eligible cells (no single treatment per cell)."
  #   )
  # }
  # keep_cells <- tmp[eligible, , drop = FALSE]
  # jdat_keep <- subset(jdat, jdat$CellID %in% keep_cells$CellID)

  # per-cell medians (y, coords)
  coord_df <- data.frame(
    sf::st_coordinates(jdat),
    sf::st_drop_geometry(jdat)
  )
  vars <- c(y, "X", "Y")
  suppressWarnings({
    med_cell <-
      coord_df |>
      dplyr::group_by(.data$.trt, .data$CellID) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(vars), ~ stats::median(.x, na.rm = TRUE)),
        .groups = "drop"
      )
  })
  med_cell <- dplyr::left_join(
    med_cell,
    unique(sf::st_drop_geometry(jdat[, c("CellID", ".trt", x)])),
    by = c("CellID", ".trt")
  )

  # model + spatial autocorr (keep spdep warnings untouched)
  fml <- stats::as.formula(paste(y, x, sep = " ~ "))
  fit <- stats::lm(fml, data = med_cell)
  med_cell$resid <- stats::residuals(fit)

  pts <- sf::st_as_sf(med_cell, coords = c("X", "Y"), crs = sf::st_crs(data))
  k1 <-
    suppressWarnings(
      spdep::knn2nb(spdep::knearneigh(sf::st_coordinates(
        pts
      )))
    )
  dmax <- max(unlist(spdep::nbdists(k1, sf::st_coordinates(pts))))
  nb <-
    suppressWarnings(
      spdep::dnearneigh(sf::st_coordinates(pts), 0, dmax)
    )
  dlist <- spdep::nbdists(nb, sf::st_coordinates(pts))
  glist <- lapply(dlist, function(x) 1 / (x / dmax))
  lw <- spdep::nb2listw(nb, glist = glist, style = "W", zero.policy = TRUE)

  MI <- spdep::moran.test(med_cell$resid, lw, zero.policy = TRUE)$estimate[[
    "Moran I statistic"
  ]]
  rho <- spatialreg::aple(med_cell$resid, lw)
  rho <- min(1, max(0, rho))

  n <- nrow(med_cell)
  ess <- round(n_eff(n = n, rho = rho), 0)

  # permutations
  trt_stats <-
    sf::st_drop_geometry(pts) |>
    dplyr::group_by(.data[[x]]) |>
    dplyr::summarise(
      !!y := stats::median(.data[[y]], na.rm = TRUE),
      .groups = "drop"
    )
  names(trt_stats)[names(trt_stats) == y] <- paste0(y, "_mean")

  pairs <- utils::combn(unique(pts$.trt), 2, simplify = TRUE)
  avail <- table(pts$.trt)
  min_avail <- min(as.integer(avail))
  min_per_group <- max(2L, min(min_avail, ceiling(ess / treatments)))

  one_run <- function(run_id) {
    base <- pts |>
      dplyr::group_by(.data$.trt) |>
      dplyr::slice_sample(n = min_per_group, replace = FALSE) |>
      dplyr::ungroup()

    per_pair <- function(j) {
      tab <- base[base$.trt %in% pairs[, j], ]
      if (any(table(tab$.trt) < 2) || stats::var(tab[[y]], na.rm = TRUE) == 0) {
        return(data.frame(
          Trt_1 = pairs[, j][1],
          Trt_2 = pairs[, j][2],
          Comparison = paste(unique(tab$.trt), collapse = " vs. "),
          p_value = NA_real_,
          run = run_id
        ))
      }
      suppressWarnings({
        a <- permuco::aovperm(fml, data = tab, np = n_p)
      })
      data.frame(
        Trt_1 = pairs[, j][1],
        Trt_2 = pairs[, j][2],
        Comparison = paste(unique(tab$.trt), collapse = " vs. "),
        p_value = a$table$`resampled P(>F)`[1],
        run = run_id
      )
    }
    do.call(rbind, lapply(seq_len(ncol(pairs)), per_pair))
  }

  perm_runs <- withr::with_seed(
    7,
    do.call(rbind, lapply(seq_len(n_s), one_run))
  )

  p_summary <-
    perm_runs |>
    dplyr::group_by(.data$Comparison) |>
    dplyr::summarise(
      p_value = stats::median(.data$p_value, na.rm = TRUE),
      Trt_1 = unique(.data$Trt_1),
      Trt_2 = unique(.data$Trt_2),
      .groups = "drop"
    ) |>
    dplyr::left_join(trt_stats, by = c("Trt_1" = x)) |>
    (function(d) d[order(d[[paste0(y, "_mean")]]), ])()

  pv <- p_summary$p_value
  names(pv) <- gsub(" vs. ", "-", p_summary$Comparison)
  pv[is.na(pv)] <- 1
  pv_adj <- stats::p.adjust(pv, method = p_adjust_method)

  letters <- multcompView::multcompLetters(pv_adj, threshold = alpha)
  letters_vec <- if ("monospacedLetters" %in% names(letters)) {
    letters$monospacedLetters
  } else {
    letters$Letters
  }
  letters_tbl <- data.frame(
    Treatment = names(letters_vec),
    .groups = unname(letters_vec),
    check.names = FALSE
  )
  names(letters_tbl)[1] <- x

  means_comp <-
    dplyr::left_join(trt_stats, letters_tbl, by = x) |>
    (function(d) d[order(d[[paste0(y, "_mean")]], decreasing = TRUE), ])()

  # report
  rep_counts <- if ("n_obs" %in% names(grid)) {
    grid$n_obs
  } else {
    lengths(sf::st_intersects(grid, data, sparse = TRUE))
  }
  gen_info <- data.frame(
    Cellsize = paste(used_params$cellsize, collapse = " x "),
    Min.Obs.Cell = suppressWarnings(min(rep_counts, na.rm = TRUE)),
    Max.Obs.Cell = suppressWarnings(max(rep_counts, na.rm = TRUE)),
    Median.Obs.Cell = suppressWarnings(stats::median(rep_counts, na.rm = TRUE)),
    Total.Cells = nrow(grid),
    # Selected.Cells = length(unique(keep_cells$CellID)),
    n = n,
    ESS = ess,
    Rho = rho,
    MoranI = MI
  )

  grid_out <- switch(
    keep_components,
    full = grid,
    light = grid[, "CellID", drop = FALSE],
    none = NULL
  )

  p_summary$p_adj <- pv_adj[gsub(" vs. ", "-", p_summary$Comparison)]

  out <- list(
    `General information` = gen_info,
    # `Cells per treatment` = table(keep_cells$trt),
    `ANOVA permutation test` = p_summary[, c("Comparison", "p_value", "p_adj")],
    `Means comparison` = means_comp,
    perm_runs = perm_runs,
    params = list(
      alpha = alpha,
      p_adjust_method = p_adjust_method,
      n_p = n_p,
      n_s = n_s,
      grid_source = grid_source,
      used_config = used_params,
      # min_per_cell_applied = min_applied,
      min_per_cell_used = min_used,
      note = info_note
    )
  )
  if (!is.null(grid_out)) {
    out$`._components` <- list(
      grid = grid_out,
      points = data[, c(x, y), drop = FALSE],
      selected_cells = dplyr::filter(grid_out, CellID %in% jdat$CellID)
    )
  }
  class(out) <- c("ofemt_result", class(out))
  out
}
