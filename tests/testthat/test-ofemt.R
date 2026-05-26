test_that("ofemt() runs end-to-end on ofe_f2 and reproduces documented numbers", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  # Suppress the internal "grid not provided" message and the spdep
  # sub-graph warnings; we test contents, not side messages.
  res <- suppressMessages(suppressWarnings(
    ofemt(
      data         = ofe_f2,
      y            = "Yield_tn",
      x            = "Treatment",
      cellsize     = 9,
      min_per_cell = 4,
      n_p          = 100,
      n_s          = 10,
      alpha        = 0.05,
      seed         = 7L
    )
  ))

  expect_s3_class(res, "ofemt_result")
  expect_true(all(c(
    "General information",
    "Cells per treatment",
    "ANOVA permutation test",
    "Means comparison",
    "perm_runs",
    "params"
  ) %in% names(res)))

  gi <- res[["General information"]]
  expect_equal(gi$Selected.Cells, 554)
  # ESS and Rho are deterministic given the data + cellsize + min_per_cell.
  expect_equal(gi$ESS, 97)
  expect_equal(gi$Rho, 0.652285, tolerance = 1e-3)

  # Moran's I — narrow tolerance since both the data and weights matrix are fixed.
  expect_equal(gi$MoranI, 0.548, tolerance = 1e-2)

  # The fertilized strip ought to come out higher in the means table.
  mc <- res[["Means comparison"]]
  expect_true("Treatment" %in% names(mc))
  expect_true("Yield_tn_mean" %in% names(mc))
  fert <- mc$Yield_tn_mean[mc$Treatment == "Fertilized"]
  ctrl <- mc$Yield_tn_mean[mc$Treatment == "Control"]
  expect_gt(fert, ctrl)
})

test_that("ofemt() is reproducible across runs with the same seed", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  args <- list(
    data         = ofe_f2,
    y            = "Yield_tn",
    x            = "Treatment",
    cellsize     = 9,
    min_per_cell = 4,
    n_p          = 100,
    n_s          = 10,
    alpha        = 0.05,
    seed         = 42L
  )

  r1 <- suppressMessages(suppressWarnings(do.call(ofemt, args)))
  r2 <- suppressMessages(suppressWarnings(do.call(ofemt, args)))

  expect_equal(
    r1[["ANOVA permutation test"]]$p_value,
    r2[["ANOVA permutation test"]]$p_value
  )
  expect_equal(r1$perm_runs$p_value, r2$perm_runs$p_value)
})

test_that("ofemt() accepts a pre-built grid and matches the auto-grid result", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  g <- make_ofe_grid(
    ofe_f2,
    x = "Treatment",
    cellsize = 9,
    min_per_cell = 4
  )

  base_args <- list(
    data  = ofe_f2,
    y     = "Yield_tn",
    x     = "Treatment",
    n_p   = 100,
    n_s   = 10,
    alpha = 0.05,
    seed  = 7L
  )

  r_auto <- suppressMessages(suppressWarnings(do.call(
    ofemt,
    c(base_args, list(cellsize = 9, min_per_cell = 4))
  )))
  r_grid <- suppressMessages(suppressWarnings(do.call(
    ofemt,
    c(base_args, list(grid = g))
  )))

  expect_equal(
    r_auto[["General information"]]$Selected.Cells,
    r_grid[["General information"]]$Selected.Cells
  )
  expect_equal(
    r_auto[["ANOVA permutation test"]]$p_value,
    r_grid[["ANOVA permutation test"]]$p_value
  )
})

test_that("ofemt() emits an informative message when grid = NULL", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  expect_message(
    suppressWarnings(
      ofemt(
        data = ofe_f2,
        y = "Yield_tn",
        x = "Treatment",
        cellsize = 9,
        min_per_cell = 4,
        n_p = 50,
        n_s = 5
      )
    ),
    "building one internally"
  )
})

test_that("ofemt() errors when the response column is non-numeric", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())
  bad <- ofe_f2
  bad$Yield_tn <- as.character(bad$Yield_tn)

  expect_error(
    ofemt(bad, y = "Yield_tn", x = "Treatment", cellsize = 9, min_per_cell = 4),
    "must be numeric"
  )
})

test_that("ofemt() errors when fewer than two treatments are present", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())
  single <- ofe_f2
  single$Treatment <- "OnlyOne"

  expect_error(
    ofemt(single, y = "Yield_tn", x = "Treatment", cellsize = 9, min_per_cell = 4),
    "two treatments"
  )
})

test_that("nmin_cell and alpha_bonferroni emit deprecation warnings", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  # Use withCallingHandlers to swallow the (unrelated) spdep sub-graph
  # warnings but let the deprecation warning bubble up to expect_warning().
  swallow_spdep <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        if (grepl("sub-graphs", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  expect_warning(
    suppressMessages(swallow_spdep(
      ofemt(
        ofe_f2,
        y = "Yield_tn",
        x = "Treatment",
        cellsize = 9,
        nmin_cell = 4,
        n_p = 50,
        n_s = 5,
        seed = 1L
      )
    )),
    "nmin_cell.*deprecated"
  )

  expect_warning(
    suppressMessages(swallow_spdep(
      ofemt(
        ofe_f2,
        y = "Yield_tn",
        x = "Treatment",
        cellsize = 9,
        min_per_cell = 4,
        n_p = 50,
        n_s = 5,
        seed = 1L,
        alpha_bonferroni = TRUE
      )
    )),
    "alpha_bonferroni.*deprecated"
  )
})
