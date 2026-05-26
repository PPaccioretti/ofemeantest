# Test helpers shared across files.
# Files prefixed with "helper-" are sourced by testthat before any test runs.

# Skip helper for tests that need the bundled ofe_f2 dataset and the heavy
# spatial dependencies. Tests that only exercise pure-R helpers (e.g. n_eff)
# should NOT use this — they should run on every platform.
skip_if_no_spatial_deps <- function() {
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("spdep")
  testthat::skip_if_not_installed("spatialreg")
  testthat::skip_if_not_installed("permuco")
  testthat::skip_if_not_installed("multcompView")
}

# Tiny synthetic OFE dataset: 2 strips of 5x10 points each on a projected CRS.
# Yields are constructed with a real treatment effect plus mild spatial
# autocorrelation so ofemt() has something signal-bearing to test against.
make_toy_ofe <- function(seed = 1L) {
  withr::with_seed(seed, {
    gx <- rep(seq(0, 90, by = 10), times = 10)
    gy <- rep(seq(0, 90, by = 10), each = 10)
    treatment <- ifelse(gx < 50, "Control", "Fertilized")
    base <- ifelse(treatment == "Fertilized", 6, 5)
    # Small AR-like noise: nearby points share a row offset.
    row_offset <- stats::rnorm(10, sd = 0.2)[as.integer(gy / 10) + 1]
    yield <- base + row_offset + stats::rnorm(length(gx), sd = 0.3)
  })

  df <- data.frame(
    X = gx,
    Y = gy,
    Treatment = treatment,
    Yield_tn = yield,
    stringsAsFactors = FALSE
  )

  # EPSG:32720 (WGS 84 / UTM 20S) — projected CRS, matches ofe_f2.
  sf::st_as_sf(df, coords = c("X", "Y"), crs = 32720)
}
