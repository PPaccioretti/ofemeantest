test_that("make_ofe_grid() returns an ofe_grid with the expected structure", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()

  g <- make_ofe_grid(pts, x = "Treatment", cellsize = 20, min_per_cell = 1)

  expect_s3_class(g, "ofe_grid")
  expect_true(all(c("grid_all", "grid_sel", "params") %in% names(g)))
  expect_s3_class(g$grid_all, "sf")
  expect_s3_class(g$grid_sel, "sf")
  expect_true("CellID" %in% names(g$grid_all))
  expect_true("CellID" %in% names(g$grid_sel))
  # Selected grid is a subset of the full grid.
  expect_true(nrow(g$grid_sel) <= nrow(g$grid_all))
})

test_that("make_ofe_grid() is deterministic for the same arguments", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()

  g1 <- make_ofe_grid(pts, x = "Treatment", cellsize = 20, min_per_cell = 1)
  g2 <- make_ofe_grid(pts, x = "Treatment", cellsize = 20, min_per_cell = 1)

  expect_identical(nrow(g1$grid_all), nrow(g2$grid_all))
  expect_identical(nrow(g1$grid_sel), nrow(g2$grid_sel))
  expect_identical(g1$grid_sel$CellID, g2$grid_sel$CellID)
})

test_that("select_grid() filters out cells with fewer than min_per_cell points", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()

  # cellsize = 8 → many small cells, most should hold 0 or 1 points.
  loose <- make_ofe_grid(pts, x = "Treatment", cellsize = 8, min_per_cell = 1)
  tight <- make_ofe_grid(pts, x = "Treatment", cellsize = 8, min_per_cell = 3)

  expect_lt(nrow(tight$grid_sel), nrow(loose$grid_sel))
  # Selected cells must all satisfy n_obs >= 3.
  expect_true(all(tight$grid_sel$n_obs >= 3))
})

test_that("select_grid() drops cells with more than one treatment", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()

  g <- make_ofe_grid(pts, x = "Treatment", cellsize = 30, min_per_cell = 1)
  # Every selected cell must contain exactly one treatment.
  expect_true(all(g$grid_sel$n_trt == 1L))
})

test_that("params carry construction arguments through to the grid object", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()

  g <- make_ofe_grid(
    pts,
    x = "Treatment",
    cellsize = 20,
    min_per_cell = 2,
    shift = c(3, 5),
    angle_deg = 0,
    buffer = 10
  )

  expect_equal(g$params$cellsize, c(20, 20))
  expect_equal(g$params$shift, c(3, 5))
  expect_equal(g$params$buffer, 10)
  expect_equal(g$params$min_per_cell, 2L)
})

test_that("filter_per_treatment() drops mixed-treatment cells", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()

  # Build a deliberately coarse grid so some cells span both treatments.
  raw <- ofemeantest:::all_cells_grid(pts, cellsize = 60)
  raw$grid_sel <- raw$grid_all # pretend nothing is filtered yet

  filtered <- ofemeantest:::filter_per_treatment(
    data = pts,
    grid = raw,
    x = "Treatment"
  )

  # All retained cells must contain points from a single treatment.
  pts_in <- sf::st_join(pts, filtered$grid_sel[, "CellID"], left = FALSE)
  per_cell <- tapply(pts_in$Treatment, pts_in$CellID, function(z) length(unique(z)))
  expect_true(all(per_cell == 1L))
})
