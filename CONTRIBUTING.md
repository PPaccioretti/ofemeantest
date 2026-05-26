# Contributing to ofemeantest

Thanks for your interest in improving the package. This document covers
the local dev setup we expect contributors to have before opening a PR.

## Development setup

1.  Clone and open the project in your R IDE:

    ``` bash
    git clone https://github.com/PPaccioretti/ofemeantest.git
    cd ofemeantest
    ```

2.  Install the dev dependencies (one-off):

    ``` r

    install.packages(c("devtools", "roxygen2", "testthat", "pkgdown"))
    ```

3.  **Enable the project’s git hooks** (one-off per clone). This points
    git at the versioned `.githooks/` directory so the pre-push check
    runs on every `git push`:

    ``` bash
    git config core.hooksPath .githooks
    ```

    You can verify it took effect with `git config --get core.hooksPath`
    (should print `.githooks`).

## What the pre-push hook does

Before every `git push` the hook runs `devtools::document()` and checks
whether anything in `man/`, `NAMESPACE` or `DESCRIPTION` changed. If it
did, the push is **blocked** with a list of files to commit. This
prevents the most common source of “stale docs” PRs — editing a roxygen
comment in `R/foo.R` and forgetting to regenerate the `.Rd`.

To bypass the hook for a single push (e.g., emergency hotfix where the
doc regen is unrelated):

``` bash
git push --no-verify
```

The same check runs in GitHub Actions on every PR
(`.github/workflows/document.yaml`), so even pushes that bypass the
local hook will be caught at review time.

## Running tests locally

``` r

devtools::test()       # run testthat suite
devtools::check()      # full R CMD check
pkgdown::build_site()  # rebuild the docs site under docs/
```

The `docs/` directory is gitignored — it’s only used for local preview
and is rebuilt by CI on push to `main`.

## Style

- Keep functions documented with roxygen2; mark internals with
  `@keywords internal`.
- Add tests for new public behaviour. Contract tests in
  `tests/testthat/test-dependency-contracts.R` guard the interfaces we
  rely on from `spdep`, `spatialreg`, `permuco` and `multcompView` — add
  a new entry there when you start depending on a new upstream function.
- Prefer base R over `dplyr`/`tidyr` to keep the Imports list small.
