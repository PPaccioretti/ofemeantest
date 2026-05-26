# Contract tests: these guard against silent breakage when upstream
# packages change the shape of objects we depend on. Each test is written
# so that a future version of spdep / spatialreg / permuco that renames or
# re-orders fields will produce a CLEAR test failure, not a silent NA.

test_that("spdep::moran.test() exposes the Moran I statistic under a stable name", {
  skip_if_no_spatial_deps()

  # Build a small spatial weights structure and a residual vector.
  pts <- make_toy_ofe()
  nb <- spdep::knn2nb(spdep::knearneigh(pts, k = 4))
  lw <- spdep::nb2listw(nb, style = "W")
  vals <- stats::rnorm(nrow(pts))

  mt <- spdep::moran.test(vals, lw)

  # ofemt() relies on this exact path:
  expect_s3_class(mt, "htest")
  expect_true(
    "estimate" %in% names(mt),
    info = paste(
      "spdep::moran.test() no longer returns $estimate.",
      "Update ofemt.R's Moran I extraction."
    )
  )
  expect_true(
    "Moran I statistic" %in% names(mt$estimate),
    info = paste(
      "spdep::moran.test()$estimate no longer contains the 'Moran I statistic'",
      "entry. The current key names are:",
      paste(shQuote(names(mt$estimate)), collapse = ", ")
    )
  )

  # And it must round-trip through the same expression we use in ofemt():
  val <- unname(mt$estimate["Moran I statistic"])
  expect_true(is.finite(val))
  expect_false(is.na(val))
})

test_that("spatialreg::aple() returns a single finite numeric", {
  skip_if_no_spatial_deps()

  pts <- make_toy_ofe()
  nb <- spdep::knn2nb(spdep::knearneigh(pts, k = 4))
  lw <- spdep::nb2listw(nb, style = "W")
  # aple() requires mean-centred input.
  resid <- stats::rnorm(nrow(pts))
  resid <- resid - mean(resid)

  rho <- spatialreg::aple(resid, lw)
  expect_length(rho, 1L)
  expect_true(is.numeric(rho))
  expect_true(is.finite(rho))
})

test_that("permuco::aovperm() returns a table with 'resampled P(>F)'", {
  skip_if_no_spatial_deps()

  # Tiny synthetic dataset so the test is fast.
  df <- data.frame(
    y = c(stats::rnorm(20, 5), stats::rnorm(20, 6)),
    g = rep(c("A", "B"), each = 20),
    stringsAsFactors = FALSE
  )
  res <- suppressWarnings(
    permuco::aovperm(y ~ g, data = df, np = 50)
  )

  expect_true(
    "table" %in% names(res),
    info = "permuco::aovperm() return changed: no $table component."
  )
  expect_true(
    "resampled P(>F)" %in% colnames(res$table),
    info = paste(
      "permuco::aovperm()$table no longer has a 'resampled P(>F)' column.",
      "Current columns:",
      paste(shQuote(colnames(res$table)), collapse = ", ")
    )
  )

  pv <- res$table[["resampled P(>F)"]][1]
  expect_true(is.finite(pv))
  expect_true(pv >= 0 && pv <= 1)
})

test_that("multcompView::multcompLetters() returns $Letters keyed by treatment", {
  pvec <- c("A-B" = 0.01, "A-C" = 0.20, "B-C" = 0.30)
  out <- multcompView::multcompLetters(pvec, threshold = 0.05)

  expect_true(
    "Letters" %in% names(out),
    info = "multcompView::multcompLetters() no longer returns $Letters."
  )
  expect_setequal(names(out$Letters), c("A", "B", "C"))
})
