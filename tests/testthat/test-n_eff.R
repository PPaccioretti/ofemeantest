# n_eff() is pure-R math; no spatial deps needed.

test_that("n_eff equals n when rho == 0", {
  expect_equal(ofemeantest:::n_eff(n = 100, rho = 0), 100)
  expect_equal(ofemeantest:::n_eff(n = 50, rho = 0), 50)
})

test_that("n_eff returns a value <= n for rho in (0, 1]", {
  n <- 200
  for (rho in c(0.1, 0.3, 0.5, 0.8, 1.0)) {
    val <- ofemeantest:::n_eff(n = n, rho = rho)
    expect_lte(val, n)
    expect_true(is.finite(val))
  }
})

test_that("n_eff is monotonically decreasing in rho", {
  n <- 100
  vals <- vapply(
    seq(0, 1, by = 0.1),
    function(r) ofemeantest:::n_eff(n = n, rho = r),
    numeric(1)
  )
  # Each step should be <= the previous one (strictly decreasing for rho > 0).
  expect_true(all(diff(vals) <= 0))
})

test_that("n_eff matches the published reference value for ofe_f2", {
  # Reference numbers from README.md / smoke test: n = 554, rho = 0.652
  # The Griffith (2005) formula should round to ~97.
  val <- ofemeantest:::n_eff(n = 554, rho = 0.652285)
  expect_equal(round(val, 0), 97)
})
