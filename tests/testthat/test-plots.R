test_that("plot_grid_selection() runs without error on an ofe_grid", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()
  g <- make_ofe_grid(pts, x = "Treatment", cellsize = 20, min_per_cell = 1)

  # Open a null PDF device so plot() side effects don't leak.
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)

  expect_silent(plot_grid_selection(g))
  expect_silent(plot_grid_selection(g, data = pts))
})

test_that("plot_pvalue_hist() returns a ggplot when ggplot2 is available", {
  skip_if_no_spatial_deps()
  testthat::skip_if_not_installed("ggplot2")

  data("ofe_f2", package = "ofemeantest", envir = environment())
  res <- suppressMessages(suppressWarnings(
    ofemt(
      ofe_f2,
      y = "Yield_tn",
      x = "Treatment",
      cellsize = 9,
      min_per_cell = 4,
      n_p = 50,
      n_s = 5,
      seed = 1L
    )
  ))

  p <- plot_pvalue_hist(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot_pvalue_hist() rejects non-ofemt_result inputs", {
  expect_error(plot_pvalue_hist(list(a = 1)), "ofemt_result")
})
