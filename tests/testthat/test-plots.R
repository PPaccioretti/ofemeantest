test_that("plot_grid_selection() runs without error on an ofe_grid", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()
  g <- make_ofe_grid(pts, x = "Treatment", cellsize = 20, min_per_cell = 1)

  # Open a null PDF device so plot() side effects don't leak.
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)

  expect_silent(plot_grid_selection(g))
  expect_silent(plot_grid_selection(g, data = pts))
  expect_silent(plot(g, data = pts))
})

test_that("plot_grid_selection() works straight off an ofemt_result", {
  skip_if_no_spatial_deps()
  pts <- make_toy_ofe()

  run <- function(keep) {
    suppressMessages(suppressWarnings(ofemt(
      pts,
      y = "Yield_tn",
      x = "Treatment",
      cellsize = 20,
      min_per_cell = 1,
      n_p = 20,
      n_s = 3,
      seed = 1L,
      keep_components = keep
    )))
  }

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)

  # "full" carries grid + points: no extra arguments needed.
  expect_silent(plot_grid_selection(run("full")))
  expect_silent(plot(run("full")))
  # "light" carries the grid only.
  expect_silent(plot_grid_selection(run("light")))
  # "none" cannot be plotted, and says so.
  expect_error(plot_grid_selection(run("none")), "keep_components")
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

test_that("plot_pvalue_hist() picks adjusted p-values when they exist", {
  skip_if_no_spatial_deps()
  testthat::skip_if_not_installed("ggplot2")

  pts <- make_toy_ofe()
  res <- suppressMessages(suppressWarnings(
    ofemt(
      pts,
      y = "Yield_tn",
      x = "Treatment",
      cellsize = 20,
      min_per_cell = 1,
      n_p = 50,
      n_s = 5,
      seed = 1L,
      p_adjust_method = "bonferroni"
    )
  ))

  p_auto <- plot_pvalue_hist(res)
  expect_s3_class(p_auto, "ggplot")
  expect_match(p_auto$labels$x, "Adjusted p-value")

  p_raw <- plot_pvalue_hist(res, which = "raw")
  expect_identical(p_raw$labels$x, "p-value")

  # Two reference lines per panel: the median and alpha.
  lines_layer <- p_auto$layers[[2]]
  expect_equal(nrow(lines_layer$data), 2L * nrow(res[["ANOVA permutation test"]]))
})

test_that("plot_pvalue_hist() rejects non-ofemt_result inputs", {
  expect_error(plot_pvalue_hist(list(a = 1)), "ofemt_result")
})
