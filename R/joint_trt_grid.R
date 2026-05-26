#' Filter OFE grid cells containing a single treatment
#'
#' Performs a spatial join between point data (`sf`) and an OFE grid,
#' identifies cells containing exactly one unique treatment, and returns
#' an updated `ofe_grid` object with only those eligible cells retained
#' in its `grid_sel` component.
#'
#' @param data An `sf` object of points containing the treatment column `x`.
#' @param grid An `ofe_grid` object (as returned by [make_grid()] or
#'   similar). Its `grid_sel` component will be used and replaced with
#'   the filtered version. Must be polygonal and share CRS with `data`.
#' @param x A character string of length one giving the name of the treatment
#'   column in `data`.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Validates the input types and column name.
#'   \item Replaces spaces in treatment labels with dots and stores them in a
#'         temporary column `".trt"`.
#'   \item Ensures that the grid has a `"CellID"` column (creates one if missing).
#'   \item Joins points to grid cells via `sf::st_intersects()`.
#'   \item Aggregates per cell and keeps only cells containing a single unique
#'         treatment (`unique_trt == 1`).
#' }
#'
#' Points outside the grid are dropped (`left = FALSE` in the join).
#'
#' @return An updated `ofe_grid` object whose `grid_sel` component contains
#'   only polygons corresponding to single-treatment cells.
#'
#' @examples
#' \dontrun{
#'   filtered_grid <- filter_per_treatment(data = pts_sf, grid = my_ofe_grid, x = "Treatment")
#'   plot(filtered_grid$grid_sel["CellID"])
#' }
#'
#' @keywords internal
filter_per_treatment <- function(data, grid, x) {
  stopifnot(inherits(grid, "ofe_grid"))
  stopifnot(inherits(data, "sf"))

  if (!is.character(x) || length(x) != 1L) {
    stop("`x` must be a character vector of length one.")
  }
  if (!(x %in% colnames(data))) {
    stop("`x` must be a valid column name in `data`.")
  }

  grid_sel <- grid$grid_sel

  data[[".trt"]] <- gsub(" ", ".", as.character(data[[x]]))
  if (!"CellID" %in% names(grid_sel)) {
    grid_sel$CellID <- seq_len(nrow(grid_sel))
  }

  # join points to cells
  jdat <- sf::st_join(data, grid_sel, join = sf::st_intersects, left = FALSE)

  tmp <- stats::aggregate(
    jdat$.trt,
    by = list(CellID = jdat$CellID),
    FUN = function(z) {
      c(
        n_trt = length(unique(z)),
        trt = ifelse(length(unique(z)) == 1, unique(z), NA),
        n = length(z)
      )
    }
  )
  tmp <- do.call(data.frame, tmp)
  names(tmp) <- c("CellID", "unique_trt", "trt", "n")

  eligible <- tmp$unique_trt == 1
  if (!any(eligible)) {
    stop("No eligible cells (no single treatment per cell).")
  }

  keep_cells <- tmp[eligible, , drop = FALSE]

  unique(keep_cells$CellID)

  grid_sel <- grid_sel[grid_sel$CellID %in% keep_cells$CellID, ]

  grid$grid_sel <- grid_sel
  grid
}
