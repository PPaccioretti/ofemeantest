#' Effective sample size under spatial autocorrelation
#'
#' @description
#' Heuristic function to adjust sample size `n` by a spatial autocorrelation parameter `rho`.
#' Based on: Griffith, D.A., & Peres-Neto, P.R. (2006). Spatial Modeling in Ecology: The Flexibility of
#' Eigenfunction Spatial Analyses. Ecology, 87(10), 2603–2613.
#'
#' @param n integer, nominal sample size.
#' @param rho numeric in [0, 1], spatial autocorrelation intensity.
#' @return numeric, effective sample size.
#' @keywords internal
n_eff <- function(n, rho) {
  a <- 1 / (1 - exp(-1.92369))
  b <- (n - 1) / n
  c <- (1 - exp(-2.12373 * rho + 0.20024 * sqrt(rho)))
  n * (1 - a * b * c)
}

#' Get permutation runs
#' @param x ofemt_result
#' @return data.frame with columns Comparison, p_value, run, etc.
#' @keywords internal
get_perm_results <- function(x) {
  stopifnot(inherits(x, "ofemt_result"))
  x$perm_runs
}

#' Evaluate an expression muffling `spdep`'s connectivity notes
#'
#' @description
#' `spdep` emits a warning whenever a neighbour object is not fully connected
#' (`"neighbour object has N sub-graphs"`). For the cell-median point pattern
#' used here that is expected and harmless: the weights matrix is built from
#' every cell's nearest-neighbour distance, so disconnected components simply
#' reflect the shape of the trial. Only the *absence* of neighbours is a real
#' problem, and that is checked explicitly by the caller (see [nb_no_links()]).
#'
#' Warnings that do not mention sub-graphs are left untouched.
#'
#' @param expr Expression to evaluate.
#' @return The value of `expr`.
#' @keywords internal
without_subgraph_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("sub-?graph", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Regions with no neighbours in an `nb` object
#'
#' In `spdep`, a region with an empty neighbour set is stored as the integer
#' `0L`. This helper returns the positions of those regions.
#'
#' @param nb An `nb` object (e.g. from `spdep::dnearneigh()`).
#' @return Integer vector of positions with an empty neighbour set.
#' @keywords internal
nb_no_links <- function(nb) {
  which(vapply(
    nb,
    function(z) length(z) == 1L && !is.na(z[1L]) && z[1L] == 0L,
    logical(1)
  ))
}

#' Quote a column name for use in a formula
#'
#' Column names coming from real field data routinely contain spaces, `+`,
#' accents or parentheses. Backticking them keeps [stats::as.formula()] happy
#' without forcing the user to rename anything.
#'
#' @param nm Character scalar, a column name.
#' @return Character scalar, the backtick-quoted name.
#' @keywords internal
quote_name <- function(nm) {
  paste0("`", gsub("`", "\\\\`", nm), "`")
}

#' Sample without the `sample()` length-one trap
#'
#' `sample(x, k)` treats a length-one numeric `x` as `1:x`. This helper always
#' samples *from the elements* of `x`.
#'
#' @param x Vector to sample from.
#' @param size Number of elements to draw.
#' @return A vector of `size` elements of `x`.
#' @keywords internal
sample_values <- function(x, size) {
  x[sample.int(length(x), size = min(size, length(x)))]
}

#' Effective sample size from distance-class Moran's I (Clifford-Dutilleul)
#'
#' @description
#' Effective sample size (ESS) for the correlation between two spatially
#' autocorrelated processes, as used by the modified *t*-test of Clifford,
#' Richardson and Hemon (1989) with the correction of Dutilleul (1993). It
#' reproduces the `ESS` component returned by `SpatialPack::modified.ttest()`
#' (up to floating-point summation order, about 1e-10 relative), re-implemented
#' here in base R so that the package does not take on a new dependency.
#'
#' @details
#' All pairwise distances between the `n` locations are split into `nclass`
#' equal-width classes spanning `[0, max distance]`. For each class \eqn{k}
#' and each process, Moran's \eqn{I_k} is computed over the pairs whose
#' distance falls in that class. The spatial correlation matrix of a process
#' is then approximated by \eqn{R_{ij} = I_{k(i,j)}} for \eqn{i \ne j} and
#' \eqn{R_{ii} = 1}, and the effective sample size is
#' \deqn{M = 1 + \frac{\mathrm{tr}(B R_x)\,\mathrm{tr}(B R_y)}{\mathrm{tr}(B R_x B R_y)},
#'   \qquad B = I - \tfrac{1}{n}\mathbf{1}\mathbf{1}'.}
#' Passing the same vector as `x` and `y` gives the ESS of a single process
#' (\eqn{R_x = R_y}), which is how [ofemt()] uses it on the ANOVA residuals.
#'
#' The \eqn{n \times n} distance matrix is never materialised: locations are
#' processed in blocks of `chunk_size` rows, so memory grows with
#' `chunk_size * n` instead of `n^2`. The result does not depend on
#' `chunk_size`. Time is \eqn{O(n^2)}, as in SpatialPack: all pairwise
#' distances have to be classified once. The maximum distance that fixes
#' the class width is taken from the convex hull, so only one pass over the
#' pairs is made. When `x` and `y` are the same vector the work for `y` is
#' skipped.
#'
#' Distance classes that contain no pair of locations get `NA` in `imoran`;
#' they do not contribute to `ESS`.
#'
#' @param x,y Numeric vectors of equal length, one value per location.
#' @param coords Numeric matrix (or data frame) with two columns of projected
#'   coordinates, one row per location.
#' @param nclass Integer, number of equal-width distance classes used to
#'   build the Moran's *I* correlogram. `SpatialPack` defaults to 13.
#' @param chunk_size Integer, number of locations processed per block. The
#'   default targets roughly 2e6 pairwise distances per block.
#' @return A list with components
#'   \describe{
#'     \item{`ESS`}{The effective sample size (numeric scalar).}
#'     \item{`upper_bounds`}{Upper bound of each distance class.}
#'     \item{`card`}{Number of location pairs in each class.}
#'     \item{`imoran`}{`nclass x 2` matrix with Moran's *I* of `x` and `y`
#'       per class.}
#'   }
#' @references
#' Clifford, P., Richardson, S. and Hemon, D. (1989). Assessing the
#' significance of the correlation between two spatial processes.
#' \emph{Biometrics}, 45, 123-134.
#'
#' Dutilleul, P. (1993). Modifying the t test for assessing the correlation
#' between two spatial processes. \emph{Biometrics}, 49, 305-314.
#'
#' Osorio, F. and Vallejos, R. (2024). SpatialPack: Tools for assessment the
#' association between two spatial processes. R package version 0.4-1.
#' \url{https://CRAN.R-project.org/package=SpatialPack}
#' @keywords internal
ess_dutilleul <- function(x, y, coords, nclass = 13L, chunk_size = NULL) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  coords <- as.matrix(coords)
  if (ncol(coords) < 2L) {
    stop("`coords` must have two columns.", call. = FALSE)
  }
  coords <- coords[, 1:2, drop = FALSE]
  storage.mode(coords) <- "double"
  n <- length(x)
  if (length(y) != n || nrow(coords) != n) {
    stop(
      "`x`, `y` and `coords` must describe the same number of locations.",
      call. = FALSE
    )
  }
  if (n < 3L) {
    stop("At least three locations are needed.", call. = FALSE)
  }
  if (anyNA(x) || anyNA(y) || anyNA(coords)) {
    stop("Missing values are not allowed.", call. = FALSE)
  }
  nclass <- as.integer(nclass)
  if (length(nclass) != 1L || is.na(nclass) || nclass < 1L) {
    stop("`nclass` must be a single positive integer.", call. = FALSE)
  }
  xc <- x - mean(x)
  yc <- y - mean(y)
  xvar <- sum(xc^2) / n
  yvar <- sum(yc^2) / n
  if (xvar <= 0 || yvar <= 0) {
    stop("`x` and `y` must not be constant.", call. = FALSE)
  }

  same_xy <- identical(x, y)

  # The largest pairwise distance fixes the class width. It is attained
  # between two vertices of the convex hull, so an O(h^2) search over the
  # hull replaces a full O(n^2) pass. `chull()` handles collinear input by
  # returning the two extreme points.
  hull <- grDevices::chull(coords)
  hx <- coords[hull, 1L]
  hy <- coords[hull, 2L]
  maxdist <- if (length(hull) < 2L) {
    0
  } else {
    max(sqrt(outer(hx, hx, "-")^2 + outer(hy, hy, "-")^2))
  }
  if (maxdist <= 0) {
    stop("All locations coincide; distance classes cannot be built.",
      call. = FALSE
    )
  }
  # Accumulated sum (not `seq_len(nclass) * width`) so that the bounds are
  # bit-for-bit those of SpatialPack.
  upper_bounds <- cumsum(rep(maxdist / nclass, nclass))

  # Single pass over the pairs, in blocks of rows. Per class we accumulate
  # the cross-product sums (for Moran's I) and, per location, the number of
  # other locations in each class, which is enough to rebuild every row sum
  # of R without storing R. Class 0 is reserved for the i == j entries.
  if (is.null(chunk_size)) {
    chunk_size <- max(1L, min(n, as.integer(2e6 %/% n)))
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  blocks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  sx <- numeric(nclass)
  sy <- numeric(nclass)
  cnt <- matrix(0, n, nclass)
  nc1 <- nclass + 1L
  cx <- coords[, 1L]
  cy <- coords[, 2L]
  for (ii in blocks) {
    m <- length(ii)
    d <- outer(cx[ii], cx, "-")^2 + outer(cy[ii], cy, "-")^2
    d <- sqrt(d)
    # First class k with d <= upper_bounds[k]; anything beyond the last
    # bound (only possible through rounding) goes to the last class.
    k <- findInterval(d, upper_bounds, left.open = TRUE) + 1L
    k[k > nclass] <- nclass
    k[(ii - 1L) * m + seq_len(m)] <- 0L # column-major index of (r, ii[r])
    if (same_xy) {
      g <- rowsum(as.vector(outer(xc[ii], xc)), k, reorder = TRUE)
    } else {
      g <- rowsum(
        cbind(as.vector(outer(xc[ii], xc)), as.vector(outer(yc[ii], yc))),
        k,
        reorder = TRUE
      )
    }
    cl <- as.integer(rownames(g))
    keep <- cl > 0L
    sx[cl[keep]] <- sx[cl[keep]] + g[keep, 1L]
    if (!same_xy) {
      sy[cl[keep]] <- sy[cl[keep]] + g[keep, 2L]
    }
    # Row r, class k -> bin (r - 1) * (nclass + 1) + k + 1; column 1 of the
    # resulting matrix holds the diagonal entries and is dropped.
    bins <- tabulate(
      (rep.int(seq_len(m), n) - 1L) * nc1 + k + 1L,
      nbins = m * nc1
    )
    cnt[ii, ] <- matrix(bins, nrow = m, byrow = TRUE)[, -1L, drop = FALSE]
  }
  if (same_xy) {
    sy <- sx
  }
  # Every unordered pair was visited twice (i, j) and (j, i).
  pairs_twice <- colSums(cnt)
  card <- pairs_twice / 2
  ix <- (sx / pairs_twice) / xvar
  iy <- (sy / pairs_twice) / yvar
  ix[pairs_twice == 0] <- NA_real_
  iy[pairs_twice == 0] <- NA_real_

  # Dutilleul (1993): M = 1 + tr(B Rx) tr(B Ry) / tr(B Rx B Ry).
  # With R_ij = I_k(i,j), R_ii = 1, the traces only need the row sums of
  # Rx and Ry and the per-class pair counts.
  ix0 <- ifelse(is.na(ix), 0, ix)
  iy0 <- ifelse(is.na(iy), 0, iy)
  rx <- 1 + as.vector(cnt %*% ix0)
  ry <- 1 + as.vector(cnt %*% iy0)
  sxx <- sum(rx)
  syy <- sum(ry)
  sxy <- sum(rx * ry)
  tr_brx <- n - sxx / n
  tr_bry <- n - syy / n
  tr_brxbry <- n + sum(pairs_twice * ix0 * iy0) +
    (sxx * syy / n - 2 * sxy) / n
  ess <- 1 + tr_brx * tr_bry / tr_brxbry

  imoran <- cbind(x = ix, y = iy)
  rownames(imoran) <- NULL
  list(
    ESS = ess,
    upper_bounds = upper_bounds,
    card = card,
    imoran = imoran
  )
}
