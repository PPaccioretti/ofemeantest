## Quick smoke test of the refactored ofemt()
## Run from the project root: Rscript data-raw/smoke_test.R

`%||%` <- function(a, b) if (is.null(a)) b else a

proj_root <- normalizePath(".")

r_files <- list.files(file.path(proj_root, "R"), pattern = "\\.R$", full.names = TRUE)
for (f in r_files) source(f, local = FALSE)

load(file.path(proj_root, "data", "ofe_f2.rda"))

cat("--- Step 1: build grid manually with make_ofe_grid() ---\n")
g <- make_ofe_grid(
  data = ofe_f2,
  x = "Treatment",
  cellsize = 9,
  min_per_cell = 4
)
cat(sprintf("Total cells: %d | Selected cells: %d\n",
            nrow(g$grid_all), nrow(g$grid_sel)))
print(g$params)

cat("\n--- Step 2: ofemt() with pre-built grid ---\n")
res1 <- ofemt(
  data       = ofe_f2,
  y          = "Yield_tn",
  x          = "Treatment",
  grid       = g,
  n_p        = 200,
  n_s        = 20,
  alpha      = 0.05
)
print(res1)

cat("\n--- Step 3: ofemt() auto-grid (grid = NULL) ---\n")
res2 <- ofemt(
  data       = ofe_f2,
  y          = "Yield_tn",
  x          = "Treatment",
  cellsize   = 9,
  min_per_cell = 4,
  n_p        = 200,
  n_s        = 20,
  alpha      = 0.05
)
print(res2)

cat("\n--- Step 4: reproducibility check (same seed -> same p-value) ---\n")
res3 <- ofemt(
  data       = ofe_f2,
  y          = "Yield_tn",
  x          = "Treatment",
  cellsize   = 9,
  min_per_cell = 4,
  n_p        = 200,
  n_s        = 20,
  alpha      = 0.05,
  seed       = 7L
)
cat("res2$p_value: ", res2$`ANOVA permutation test`$p_value, "\n")
cat("res3$p_value: ", res3$`ANOVA permutation test`$p_value, "\n")
cat("Equal?       ", isTRUE(all.equal(
  res2$`ANOVA permutation test`$p_value,
  res3$`ANOVA permutation test`$p_value
)), "\n")

cat("\n--- Step 5: class check ---\n")
cat("class(res1): ", paste(class(res1), collapse = ", "), "\n")
cat("inherits ofemt_result? ", inherits(res1, "ofemt_result"), "\n")

cat("\n--- Step 6: perm_runs structure (for plot_pvalue_hist) ---\n")
str(res1$perm_runs)

cat("\nSmoke test completed successfully.\n")
