# Building a ggplot is lazy: errors in the layers, scales or guides only
# surface when the plot is rendered. Every ggplot assertion below therefore
# goes through ggplot_build(), otherwise a broken plot would pass silently.
expect_builds <- function(p) {
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_no_error(ggplot2::ggplot_build(p))
}

test_that("plot_grid_selection() runs without error on an ofe_grid", {
  skip_if_no_spatial_deps()
  testthat::skip_if_not_installed("ggplot2")
  pts <- make_toy_ofe()
  g <- make_ofe_grid(pts, x = "Treatment", cellsize = 20, min_per_cell = 1)

  expect_builds(plot_grid_selection(g, data = pts))
  expect_builds(plot(g, data = pts))
  expect_builds(plot_grid_selection(g, points = FALSE))
  # Built without return_points and given no data: say so instead of
  # silently drawing a grid with no observations on it.
  expect_message(plot_grid_selection(g), "no observations to draw")
})

test_that("an ofe_grid built with return_points = TRUE draws its points", {
  skip_if_no_spatial_deps()
  testthat::skip_if_not_installed("ggplot2")
  pts <- make_toy_ofe()

  g_pts <- make_ofe_grid(
    pts,
    x = "Treatment",
    cellsize = 20,
    min_per_cell = 1,
    return_points = TRUE
  )
  expect_false(is.null(g_pts$points_sel))

  # points_sel is picked up automatically -> no "no observations" message
  expect_silent(p <- plot_grid_selection(g_pts))
  expect_builds(p)
  expect_builds(plot(g_pts))
  # ... and points = FALSE still suppresses the layer without complaining
  expect_silent(p_nopts <- plot_grid_selection(g_pts, points = FALSE))
  expect_builds(p_nopts)
  # One layer fewer without the points
  expect_equal(length(p$layers), length(p_nopts$layers) + 1L)
})

test_that("plot_grid_selection() works straight off an ofemt_result", {
  skip_if_no_spatial_deps()
  testthat::skip_if_not_installed("ggplot2")
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

  # "full" carries grid + points: no extra arguments needed.
  expect_builds(plot_grid_selection(run("full")))
  expect_builds(plot(run("full")))
  # "light" carries the grid only, and says so when points are requested.
  expect_message(plot_grid_selection(run("light")), "no observations to draw")
  expect_builds(plot_grid_selection(run("light"), points = FALSE))
  # "none" cannot be plotted, and says so.
  expect_error(plot_grid_selection(run("none")), "keep_components")
})

test_that("the engine argument switches between ggplot2 and base", {
  skip_if_no_spatial_deps()
  testthat::skip_if_not_installed("ggplot2")
  pts <- make_toy_ofe()
  g <- make_ofe_grid(
    pts,
    x = "Treatment",
    cellsize = 20,
    min_per_cell = 1,
    return_points = TRUE
  )

  # ggplot2 is the default and returns an object
  expect_builds(plot_grid_selection(g))

  # base draws and returns NULL invisibly
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_null(plot_grid_selection(g, engine = "base"))
  expect_null(plot(g, engine = "base"))

  expect_error(plot_grid_selection(g, engine = "nope"), "'arg'")
})

test_that("legend placement works for both engines", {
  skip_if_no_spatial_deps()
  testthat::skip_if_not_installed("ggplot2")
  pts <- make_toy_ofe()
  g <- make_ofe_grid(
    pts,
    x = "Treatment",
    cellsize = 20,
    min_per_cell = 1,
    return_points = TRUE
  )

  # Base-style corner keywords are translated to ggplot2 sides, so the same
  # legend_pos value is valid whichever engine is in use.
  expect_identical(legend_pos_gg(NULL), "right")
  expect_identical(legend_pos_gg("topleft"), "left")
  expect_identical(legend_pos_gg("bottomright"), "right")
  expect_identical(legend_pos_gg("bottom"), "bottom")
  expect_identical(legend_pos_gg("none"), "none")
  expect_identical(legend_pos_gg("topleft", legend = FALSE), "none")
  expect_identical(legend_pos_gg(c(0.1, 0.9)), c(0.1, 0.9))

  expect_builds(plot_grid_selection(g, legend_pos = "bottom"))
  expect_builds(plot_grid_selection(g, legend = FALSE))

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_null(plot_grid_selection(g, engine = "base", legend = FALSE))
  expect_null(plot_grid_selection(g, engine = "base", legend_pos = "bottomright"))
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

  expect_builds(plot_pvalue_hist(res))

  # base engine draws instead of returning an object
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_null(plot_pvalue_hist(res, engine = "base"))
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
  expect_builds(p_auto)
  expect_match(p_auto$labels$x, "Adjusted p-value")

  p_raw <- plot_pvalue_hist(res, which = "raw")
  expect_identical(p_raw$labels$x, "p-value")

  # Two reference lines per panel: the median and alpha.
  lines_layer <- p_auto$layers[[2]]
  expect_equal(
    nrow(lines_layer$data),
    2L * nrow(res[["ANOVA permutation test"]])
  )
})

test_that("plot_pvalue_hist() rejects non-ofemt_result inputs", {
  expect_error(plot_pvalue_hist(list(a = 1)), "ofemt_result")
})
