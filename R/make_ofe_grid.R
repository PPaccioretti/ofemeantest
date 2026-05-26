#' Create an OFE grid and select valid cells
#'
#' Convenience wrapper that builds a grid with [all_cells_grid()] and
#' filters cells using [select_grid()]. It is the main entry point for
#' generating grids for on-farm experiment (OFE) analysis.
#'
#' @param data An `sf` object of points. See [all_cells_grid()].
#' @param x Character scalar with the treatment column name in `data`.
#'   See [select_grid()].
#' @param cellsize Grid cell size. See [all_cells_grid()].
#' @param min_per_cell Minimum observations per cell. See [select_grid()].
#' @param angle_deg Grid rotation in degrees. See [all_cells_grid()].
#' @param buffer Buffer applied to the data bounding box before gridding.
#'   See [all_cells_grid()].
#' @param shift Offset applied to the grid origin. See [all_cells_grid()].
#' @param return_points Logical. See [select_grid()].
#'
#' @return An `ofe_grid` object returned by [select_grid()], containing
#'   the full grid, the selected cells, and associated metadata.
#'
#' @seealso [all_cells_grid()], [select_grid()]
#'
#' @export
make_ofe_grid <- function(
  data,
  x,
  cellsize,
  min_per_cell = 1L,
  angle_deg = 0,
  buffer = 0,
  shift = c(0, 0),
  return_points = FALSE
) {
  the_grid <- all_cells_grid(
    data = data,
    cellsize = cellsize,
    angle_deg = angle_deg,
    buffer = buffer,
    shift = shift
  )

  the_selected_grid <- select_grid(
    grid = the_grid,
    data = data,
    x = x,
    min_per_cell = min_per_cell,
    return_points = return_points
  )
  class(the_selected_grid) <- c("ofe_grid", "ofe_grid_filtered", "list")

  the_selected_grid
}


#' Build OFE grid and (initially) select all cells
#'
#' @param data An `sf` object of points with geometry.
#' @param cellsize numeric length-1 or length-2, grid cell size.
#' @param angle_deg numeric, grid rotation in degrees.
#' @param buffer numeric, buffer added around data bbox before gridding (same units as CRS).
#' @param shift numeric length-2, offset added to grid origin (xmin, ymin).
#'
#' @return An object of class `ofe_grid` (a list) with:
#'   \itemize{
#'     \item `grid_all`: `sf` polygons of the full grid (with `CellID`).
#'     \item `grid_sel`: initially identical to `grid_all` (no filtering yet).
#'     \item `params`: a named list of grid parameters (cellsize, shift, angle, buffer, CRS).
#'   }
#' @keywords internal
all_cells_grid <- function(
  data,
  cellsize,
  angle_deg = 0,
  buffer = 0,
  shift = c(0, 0)
) {
  stopifnot(inherits(data, "sf"))
  crs <- sf::st_crs(data)
  if (isTRUE(sf::st_is_longlat(crs))) {
    stop("Use a projected CRS (not lon/lat).")
  }

  cs <- if (length(cellsize) == 1) c(cellsize, cellsize) else cellsize[1:2]
  shift <- if (length(shift) == 1) c(shift, shift) else shift[1:2]

  bb <- sf::st_bbox(data)
  if (buffer > 0) {
    bb <- c(
      bb["xmin"] - buffer,
      bb["ymin"] - buffer,
      bb["xmax"] + buffer,
      bb["ymax"] + buffer
    ) |>
      sf::st_bbox(crs = crs)
  }

  offset <- c(bb["xmin"] + shift[1], bb["ymin"] + shift[2])

  grid <- sf::st_make_grid(
    x = sf::st_as_sfc(bb, crs = crs),
    cellsize = cs,
    offset = offset,
    what = "polygons",
    square = TRUE
  ) |>
    sf::st_as_sf()
  sf::st_crs(grid) <- crs

  if (!isTRUE(all.equal(angle_deg, 0))) {
    grid <- .rotate_sfc(grid, angle_deg)
  }

  grid$CellID <- seq_len(nrow(grid))

  # compute observation counts for reference (using provided data extent)
  idx <- sf::st_intersects(grid, data, sparse = TRUE)
  grid$n_obs <- lengths(idx)

  params <- list(
    cellsize = cs,
    shift = shift,
    angle_deg = angle_deg,
    buffer = buffer,
    crs_wkt = sf::st_crs(data)$wkt %||% as.character(sf::st_crs(data))
  )

  out <- list(
    grid_all = grid,
    grid_sel = grid, # no selection yet
    params = params
  )
  class(out) <- c("ofe_grid", "list")
  out
}

#' Filter grid by minimum observations and minimum number of treatments
#'
#' Spatially joins points to the grid, counts observations and distinct treatments
#' per cell, and keeps only cells that satisfy both thresholds.
#'
#' @param grid An `ofe_grid` (from `make_grid()`) or an `sf` polygons grid with `CellID`.
#' @param data An `sf` object of points.
#' @param x Character scalar with the name of the treatment column in `data`.
#' @param min_per_cell Integer, minimum number of observations per cell (default 1).
#' @param return_points Logical, if `TRUE` also returns the joined points restricted
#'   to the selected cells (component `points_sel`). Default `FALSE`.
#'
#' @return An `ofe_grid` (list) with updated `grid_sel` and `params`. Adds:
#'   \itemize{
#'     \item `cell_stats`: a data.frame with per-cell `n_obs`, `n_trt`, and `treatments` list.
#'     \item `points_sel` (optional): the points that fall in the selected cells.
#'   }
#' @keywords internal
select_grid <- function(
  grid,
  data,
  x,
  min_per_cell = 1L,
  return_points = FALSE
) {
  stopifnot(inherits(data, "sf"), is.character(x), length(x) == 1)

  # Unpack grid (accept ofe_grid or plain sf)
  if (inherits(grid, "ofe_grid")) {
    grid_sf <- grid$grid_all
    params <- grid$params
  } else {
    grid_sf <- grid
    params <- list()
  }
  stopifnot(inherits(grid_sf, "sf"), "CellID" %in% names(grid_sf))

  # Normalize treatment labels and ensure column exists
  if (!x %in% names(data)) {
    stop(sprintf("Column '%s' not found in data.", x))
  }
  data$.trt <- gsub("\\s+", ".", data[[x]])

  # CRS check
  if (
    !identical(sf::st_crs(grid_sf)$epsg, sf::st_crs(data)$epsg) &&
      !identical(sf::st_crs(grid_sf)$wkt, sf::st_crs(data)$wkt)
  ) {
    stop("CRS mismatch between grid and data.")
  }

  # Join: bring CellID to points
  pts_join <- sf::st_join(
    data,
    grid_sf[, "CellID"],
    left = FALSE # keep only points that fall in some cell
  )

  # Per-cell stats: split treatment vector by CellID
  trt_by_cell <- split(
    sf::st_drop_geometry(pts_join)$.trt,
    sf::st_drop_geometry(pts_join)$CellID
  )
  cell_stats <- data.frame(
    CellID = as.integer(names(trt_by_cell)),
    n_obs  = lengths(trt_by_cell),
    n_trt  = vapply(trt_by_cell, function(z) length(unique(z)), integer(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  cell_stats$treatments <- unname(lapply(trt_by_cell, function(z) sort(unique(z))))

  # Merge counts back to grid (preserve grid_all row order via match())
  grid_all <- grid_sf
  grid_all$n_obs <- NULL
  m <- match(grid_all$CellID, cell_stats$CellID)
  grid_all$n_obs      <- cell_stats$n_obs[m]
  grid_all$n_trt      <- cell_stats$n_trt[m]
  grid_all$treatments <- cell_stats$treatments[m]
  grid_all$n_obs[is.na(grid_all$n_obs)] <- 0L
  grid_all$n_trt[is.na(grid_all$n_trt)] <- 0L

  keep_ids <- cell_stats$CellID[
    cell_stats$n_obs >= min_per_cell & cell_stats$n_trt == 1
  ]

  grid_sel <- grid_all[grid_all$CellID %in% keep_ids, , drop = FALSE]

  # Update params
  params <- utils::modifyList(
    params,
    list(
      min_per_cell = as.integer(min_per_cell),
      filter_note = "Selected cells have >= min_per_cell obs AND 1 distinct treatments"
    )
  )

  out <- list(
    grid_all = grid_all,
    grid_sel = grid_sel,
    params = params,
    cell_stats = cell_stats
  )

  if (isTRUE(return_points)) {
    out$points_sel <- pts_join[pts_join$CellID %in% keep_ids, , drop = FALSE]
  }

  class(out) <- c("ofe_grid", "ofe_grid_filtered", "list")
  out
}

# ------- helpers -------

.rotate_sfc <- function(x, angle_deg = 0, center = NULL) {
  if (abs(angle_deg) < 1e-12) {
    return(x)
  }
  crs_x <- sf::st_crs(x)
  bb <- sf::st_bbox(x)
  if (is.null(center)) {
    center <- c((bb["xmin"] + bb["xmax"]) / 2, (bb["ymin"] + bb["ymax"]) / 2)
  }
  th <- angle_deg * pi / 180
  R <- matrix(c(cos(th), -sin(th), sin(th), cos(th)), 2, byrow = TRUE)
  g <- sf::st_geometry(x)
  g <- (g - center) * R + center
  g <- sf::st_set_crs(g, crs_x)
  sf::st_set_geometry(x, g)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
