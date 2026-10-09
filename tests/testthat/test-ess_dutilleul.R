# ess_dutilleul() is pure base R; the SpatialPack comparisons are skipped
# when that package is not installed.

# Small irregular point set with two mildly autocorrelated variables.
make_ess_fixture <- function(n = 60L, seed = 3L) {
  withr::with_seed(seed, {
    coords <- cbind(stats::runif(n, 0, 100), stats::runif(n, 0, 100))
    trend <- coords[, 1] / 100
    x <- trend + stats::rnorm(n, sd = 0.5)
    y <- -trend + stats::rnorm(n, sd = 0.5)
  })
  list(x = x, y = y, coords = coords)
}

# Reference implementation straight from Dutilleul (1993): build the full
# n x n correlation matrices from the correlogram and take the traces.
ess_by_traces <- function(x, y, coords, nclass) {
  n <- length(x)
  d <- as.matrix(stats::dist(coords))
  ub <- cumsum(rep(max(d) / nclass, nclass))
  k <- matrix(pmin(findInterval(d, ub, left.open = TRUE) + 1L, nclass), n)
  xc <- x - mean(x)
  yc <- y - mean(y)
  moran <- function(v) {
    vapply(seq_len(nclass), function(cl) {
      sel <- which(k == cl & upper.tri(d), arr.ind = TRUE)
      if (nrow(sel) == 0) return(0)
      mean(v[sel[, 1]] * v[sel[, 2]]) / (sum(v^2) / n)
    }, numeric(1))
  }
  build_r <- function(I) {
    r <- matrix(I[k], n)
    diag(r) <- 1
    r
  }
  rx <- build_r(moran(xc))
  ry <- build_r(moran(yc))
  b <- diag(n) - matrix(1 / n, n, n)
  tr <- function(m) sum(diag(m))
  1 + tr(b %*% rx) * tr(b %*% ry) / tr(b %*% rx %*% b %*% ry)
}

test_that("ess_dutilleul() reproduces the Dutilleul (1993) trace formula", {
  f <- make_ess_fixture()
  for (nclass in c(5L, 13L)) {
    got <- ofemeantest:::ess_dutilleul(f$x, f$y, f$coords, nclass = nclass)
    expect_equal(
      got$ESS,
      ess_by_traces(f$x, f$y, f$coords, nclass),
      tolerance = 1e-10
    )
    expect_length(got$upper_bounds, nclass)
    expect_length(got$card, nclass)
    expect_equal(dim(got$imoran), c(nclass, 2L))
    # Every unordered pair lands in exactly one class.
    expect_equal(sum(got$card), length(f$x) * (length(f$x) - 1) / 2)
  }
})

test_that("ess_dutilleul() does not depend on the block size", {
  f <- make_ess_fixture()
  whole <- ofemeantest:::ess_dutilleul(f$x, f$y, f$coords)
  for (chunk in c(1L, 7L, 25L, 1000L)) {
    blocked <- ofemeantest:::ess_dutilleul(
      f$x, f$y, f$coords, chunk_size = chunk
    )
    expect_equal(blocked, whole)
  }
})

test_that("ess_dutilleul() equals n when the correlogram carries no information", {
  # Alternating +1/-1 on a line with a single distance class: Moran's I of
  # the only class is -1/(n - 1), so R = (1 + 1/(n - 1)) I - J/(n - 1), and
  # B R = (1 + 1/(n - 1)) B. The centred traces then give ESS = n exactly.
  n <- 8L
  coords <- cbind(seq_len(n), 0)
  x <- rep(c(1, -1), times = n / 2)
  got <- ofemeantest:::ess_dutilleul(x, x, coords, nclass = 1L)
  expect_equal(got$ESS, n, tolerance = 1e-12)
  expect_equal(unname(got$imoran[1, 1]), -1 / (n - 1))
})

test_that("ess_dutilleul() matches SpatialPack::modified.ttest()", {
  skip_if_not_installed("SpatialPack")

  f <- make_ess_fixture(n = 120L, seed = 11L)
  for (nclass in c(6L, 13L, 20L)) {
    ref <- SpatialPack::modified.ttest(f$x, f$y, f$coords, nclass = nclass)
    got <- ofemeantest:::ess_dutilleul(f$x, f$y, f$coords, nclass = nclass)
    expect_equal(got$ESS, ref$ESS, tolerance = 1e-10)
    expect_equal(got$upper_bounds, ref$upper.bounds, tolerance = 1e-12)
    expect_equal(got$card, ref$card)
    expect_equal(unname(got$imoran[, 1]), ref$imoran[[1]], tolerance = 1e-10)
    expect_equal(unname(got$imoran[, 2]), ref$imoran[[2]], tolerance = 1e-10)
  }

  # Same process twice, which is how ofemt() uses it.
  ref <- SpatialPack::modified.ttest(f$x, f$x, f$coords)
  got <- ofemeantest:::ess_dutilleul(f$x, f$x, f$coords)
  expect_equal(got$ESS, ref$ESS, tolerance = 1e-10)
})

test_that("ess_dutilleul() matches SpatialPack on its murray example", {
  skip_if_not_installed("SpatialPack")
  murray <- NULL
  utils::data("murray", package = "SpatialPack", envir = environment())
  skip_if(is.null(murray))

  ref <- SpatialPack::modified.ttest(
    murray$As, murray$Pb, murray[c("xpos", "ypos")]
  )
  got <- ofemeantest:::ess_dutilleul(
    murray$As, murray$Pb, murray[c("xpos", "ypos")]
  )
  expect_equal(got$ESS, ref$ESS, tolerance = 1e-10)
  expect_equal(got$card, ref$card)
})

test_that("ess_dutilleul() validates its inputs", {
  f <- make_ess_fixture(n = 10L)
  expect_error(
    ofemeantest:::ess_dutilleul(f$x, f$y[-1], f$coords),
    "same number of locations"
  )
  expect_error(
    ofemeantest:::ess_dutilleul(f$x, f$y, f$coords[, 1, drop = FALSE]),
    "two columns"
  )
  expect_error(
    ofemeantest:::ess_dutilleul(f$x, f$y, f$coords, nclass = 0),
    "nclass"
  )
  expect_error(
    ofemeantest:::ess_dutilleul(rep(1, 10), f$y, f$coords),
    "constant"
  )
  xna <- f$x
  xna[1] <- NA
  expect_error(
    ofemeantest:::ess_dutilleul(xna, f$y, f$coords),
    "Missing values"
  )
})
