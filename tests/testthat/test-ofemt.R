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

test_that("nmin_cell and alpha_bonferroni are gone", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  expect_error(
    ofemt(ofe_f2, y = "Yield_tn", x = "Treatment", cellsize = 9, nmin_cell = 4),
    "unused argument"
  )
  expect_error(
    ofemt(
      ofe_f2,
      y = "Yield_tn",
      x = "Treatment",
      cellsize = 9,
      alpha_bonferroni = TRUE
    ),
    "unused argument"
  )
})

test_that("ofemt() does not leak spdep's sub-graph connectivity warnings", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  warnings_seen <- character()
  withCallingHandlers(
    suppressMessages(ofemt(
      data = ofe_f2,
      y = "Yield_tn",
      x = "Treatment",
      cellsize = 9,
      min_per_cell = 4,
      n_p = 20,
      n_s = 3,
      seed = 7L
    )),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_false(any(grepl("sub-?graph", warnings_seen)))
})

test_that("treatment labels with special characters survive intact", {
  skip_if_no_spatial_deps()

  pts <- make_toy_ofe()
  # Labels with a space, a "+", a "-" and parentheses: every separator
  # multcompView and formula parsing could trip on.
  pts$Treatment <- ifelse(
    pts$Treatment == "Control",
    "Testigo (0 kg)",
    "Bortrac + Zin-trac"
  )

  res <- suppressMessages(suppressWarnings(ofemt(
    pts,
    y = "Yield_tn",
    x = "Treatment",
    cellsize = 20,
    min_per_cell = 1,
    n_p = 50,
    n_s = 5,
    seed = 1L
  )))

  mc <- res[["Means comparison"]]
  expect_setequal(mc$Treatment, c("Testigo (0 kg)", "Bortrac + Zin-trac"))
  expect_false(anyNA(mc$letters))
  expect_equal(
    res[["ANOVA permutation test"]]$Comparison,
    "Bortrac + Zin-trac vs. Testigo (0 kg)"
  )
  expect_setequal(names(res[["Cells per treatment"]]), mc$Treatment)
})

test_that("compact letters follow the ordering of the means", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  res <- suppressMessages(suppressWarnings(ofemt(
    data = ofe_f2,
    y = "Yield_tn",
    x = "Treatment",
    cellsize = 9,
    min_per_cell = 4,
    n_p = 100,
    n_s = 10,
    seed = 7L
  )))

  mc <- res[["Means comparison"]]
  # Means come out in decreasing order ...
  expect_false(is.unsorted(rev(mc$Yield_tn_mean)))
  # ... and the first (highest) group always carries the first letter.
  expect_true(grepl("a", mc$letters[1], fixed = TRUE))
})

test_that("p_adj is the median of the per-run adjusted p-values", {
  skip_if_no_spatial_deps()
  data("ofe_f2", package = "ofemeantest", envir = environment())

  res <- suppressMessages(suppressWarnings(ofemt(
    data = ofe_f2,
    y = "Yield_tn",
    x = "Treatment",
    cellsize = 9,
    min_per_cell = 4,
    n_p = 100,
    n_s = 10,
    seed = 7L,
    p_adjust_method = "bonferroni"
  )))

  expect_true("p_adj" %in% names(res$perm_runs))
  tbl <- res[["ANOVA permutation test"]]
  by_hand <- vapply(
    tbl$Comparison,
    function(cmp) {
      stats::median(res$perm_runs$p_adj[res$perm_runs$Comparison == cmp])
    },
    numeric(1)
  )
  expect_equal(tbl$p_adj, unname(by_hand))
})

test_that("keep_components controls the embedded geometries", {
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

  none <- run("none")
  expect_null(none$grid)
  expect_null(none$cell_medians)
  expect_null(none$points_joined)

  light <- run("light")
  expect_s3_class(light$grid, "ofe_grid")
  expect_s3_class(light$grid$grid_all, "sf")
  expect_s3_class(light$grid$grid_sel, "sf")
  expect_null(light$cell_medians)
  expect_null(light$points_joined)

  full <- run("full")
  expect_s3_class(full$grid, "ofe_grid")
  expect_s3_class(full$cell_medians, "sf")
  expect_s3_class(full$points_joined, "sf")
  expect_true("residuos" %in% names(full$cell_medians))
  expect_true("CellID" %in% names(full$points_joined))
})
