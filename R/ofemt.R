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
#' @param grid optional: an `ofe_grid` (e.g., from [make_ofe_grid()]). If `NULL`,
#'   the grid is generated internally using `cellsize`, `min_per_cell`,
#'   `shift`, `angle_deg` and `buffer`, and an informative message is emitted.
#' @param crs optional target projected CRS if `data` is in lon/lat.
#' @param seed integer; controls reproducibility of the permutation sampling
#'   (the grid itself is deterministic from its construction arguments). Set
#'   `seed = NULL` to let results vary across runs.
#' @param nmin_cell,alpha_bonferroni **Deprecated.**
#'   Use `min_per_cell` instead of `nmin_cell`, and `p_adjust_method = "bonferroni"`
#'   instead of `alpha_bonferroni`.
#'
#' @return An object of class `ofemt_result`.
#' @export
#'
#' @details
#' The methodology involves:
#' \enumerate{
#'   \item calculation of effective sample size (ESS) given the underlying spatial
#' structure.
#'   \item ANOVA permutation test on a random sample of ESS.
#'   \item generation of the empirical distribution of p-values from repetition of
#' step two. The median of this empirical distribution is regarded as the
#' p-value associated with the non-treatment effect hypothesis.
#' }
#' The test can be easily extended to cover scenarios with more than two
#' treatments by employing pairwise comparisons of OFE across multiple
#' treatments, with p-values can be adjusted for multiplicity
#' using Bonferroni correction
#'
#' @return a list of length 4 with mean test comparisson results
#'
#' @references A new method to compare
#' treatments in unreplicated on-farm experimentation. Córdoba M.,
#' Paccioretti P., Balzarini M. Under review.
#'
#' @importFrom magrittr %>%
#' @export
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
  seed = 7L,
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
    message(
      "`grid` not provided: building one internally via `make_ofe_grid()`. ",
      "Pass a pre-built `ofe_grid` to inspect or reuse the selection."
    )
    gint <- make_ofe_grid(
      data = data,
      x = x,
      cellsize = used_params$cellsize,
      min_per_cell = used_params$min_per_cell,
      angle_deg = used_params$angle_deg,
      buffer = used_params$buffer,
      shift = used_params$shift
    )
    grid <- gint$grid_sel
    used_params <- gint$params
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

  # Cell medians (per-cell median yield + coordinates, keep treatment label)
  coords <- sf::st_coordinates(jdat)
  jdat_df <- data.frame(
    sf::st_drop_geometry(jdat),
    X = coords[, "X"],
    Y = coords[, "Y"],
    stringsAsFactors = FALSE
  )

  cell_keys <- list(CellID = jdat_df$CellID, .trt = jdat_df$.trt)
  cell_med <- stats::aggregate(
    jdat_df[, c(y, "X", "Y")],
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
  trt_map <- unique(jdat_df[, c(".trt", x)])
  cell_med <- merge(cell_med, trt_map, by = ".trt", all.x = TRUE)

  # One-way ANOVA on cell medians to get residuals for spatial diagnostics
  my_model <- stats::lm(
    stats::as.formula(paste(y, "~", ".trt")),
    data = cell_med
  )
  cell_med$residuos <- stats::residuals(my_model)

  cell_sf <- sf::st_as_sf(
    cell_med,
    coords = c("X", "Y"),
    crs = sf::st_crs(data)
  )

  # Spatial weights: nearest neighbours, distance-weighted
  k1 <- spdep::knn2nb(spdep::knearneigh(cell_sf))
  dmax <- max(unlist(spdep::nbdists(k1, cell_sf)))
  gri <- spdep::dnearneigh(cell_sf, 0, dmax)
  dist <- spdep::nbdists(gri, cell_sf)
  fdist <- lapply(dist, function(d) 1 / (d / dmax))
  lw <- tryCatch(
    spdep::nb2listw(gri, glist = fdist, style = "W"),
    error = function(e) {
      spdep::set.ZeroPolicyOption(TRUE)
      spdep::nb2listw(gri, glist = fdist, style = "W")
    }
  )

  # Spatial autocorrelation of residuals + effective sample size
  rho <- spatialreg::aple(cell_sf$residuos, lw)
  rho <- min(1, max(0, rho))
  moran_test <- spdep::moran.test(cell_sf$residuos, lw)
  MI <- unname(moran_test$estimate["Moran I statistic"])
  n <- nrow(cell_sf)
  ess <- round(n_eff(n = n, rho = rho), 0)

  # Treatment medians (on cell medians) — keep both labels for joining
  trt_med_df <- stats::aggregate(
    cell_med[[y]],
    by = list(.trt = cell_med$.trt),
    FUN = stats::median
  )
  names(trt_med_df) <- c(".trt", y)

  trt_levels <- sort(unique(cell_med$.trt))
  n_trt <- length(trt_levels)
  compar <- utils::combn(trt_levels, 2, simplify = TRUE)

  # Multiple-permutation runs
  multipermutacion <- function(p) {
    # Sample ceiling(ess / n_trt) cells per treatment without replacement
    base_perm_idx <- unlist(lapply(trt_levels, function(t) {
      idx <- which(cell_sf$.trt == t)
      sample(idx, size = min(length(idx), ceiling(ess / n_trt)))
    }))
    base_perm <- cell_sf[base_perm_idx, ]

    permt_trat <- function(j) {
      pair <- compar[, j]
      sub <- base_perm[base_perm$.trt %in% pair, ]
      suppressWarnings(
        perm_res <- permuco::aovperm(
          stats::as.formula(paste(y, "~ .trt")),
          data = sub,
          np = n_p
        )
      )
      data.frame(
        Trt_1 = pair[1],
        Trt_2 = pair[2],
        Comparison = paste(pair, collapse = " vs. "),
        p_value = perm_res$table$`resampled P(>F)`[1],
        stringsAsFactors = FALSE
      )
    }

    out <- do.call(rbind, lapply(seq_len(ncol(compar)), permt_trat))
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

  # Median p-value per comparison + multiplicity adjustment
  pvals_by_comp <- stats::aggregate(
    perm_runs$p_value,
    by = list(Comparison = perm_runs$Comparison),
    FUN = function(z) stats::median(z, na.rm = TRUE)
  )
  names(pvals_by_comp) <- c("Comparison", "p_value")
  pvals_by_comp$p_adj <- stats::p.adjust(
    pvals_by_comp$p_value,
    method = p_adjust_method
  )

  # Compact letter display from adjusted p-values
  pvec <- pvals_by_comp$p_adj
  names(pvec) <- gsub(" vs. ", "-", pvals_by_comp$Comparison)
  letras_comp <- multcompView::multcompLetters(pvec, threshold = alpha)
  my_letters <- if (!is.null(letras_comp$monospacedLetters)) {
    letras_comp$monospacedLetters
  } else {
    letras_comp$Letters
  }
  letras <- data.frame(
    .trt = names(letras_comp$Letters),
    letters = my_letters,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  means_table <- merge(trt_med_df, letras, by = ".trt")
  means_table <- means_table[order(means_table[[y]], decreasing = TRUE), ]
  # Rename to the user-facing column names
  names(means_table)[names(means_table) == ".trt"] <- x
  names(means_table)[names(means_table) == y] <- paste0(y, "_mean")
  row.names(means_table) <- NULL

  cells_per_trt <- table(cell_med[[x]])

  # Use the selected-cell n_obs as the per-cell observation counts
  n_obs_sel <- cell_med$n_obs
  infogral <- data.frame(
    Cellsize = paste(used_params$cellsize, collapse = "_"),
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
      grid_source = grid_source
    )
  )

  if (keep_components != "none") {
    result$grid <- grid
    if (keep_components == "full") {
      result$cell_medians <- cell_sf
      result$points_joined <- jdat
    }
  }

  class(result) <- c("ofemt_result", "list")
  result
}
