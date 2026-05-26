# Methodology

This vignette describes the statistical protocol implemented by the
`ofemeantest` package, intended for readers who want to understand
*what* the function computes — beyond the user-facing recipe in
[`vignette("getting-started")`](https://ppaccioretti.github.io/ofemeantest/articles/getting-started.md).

## 1. Problem setting

On-farm experiments (OFE) on commercial fields are typically laid out as
side-by-side strips with **no replication** of treatments and no
randomisation. Yield-monitor data over these strips are dense, irregular
and **strongly spatially autocorrelated**, so classical ANOVA on the raw
observations grossly inflates the effective sample size and breaks its
independence assumption.

The protocol implemented here adapts ANOVA to this setting by (i)
aggregating the raw observations into a regular grid, (ii) estimating
the autocorrelation-corrected effective sample size, and (iii) inferring
treatment effects from a distribution of permutational *p*-values
computed on random subsamples of that ESS.

## 2. Grid construction

[`make_ofe_grid()`](https://ppaccioretti.github.io/ofemeantest/reference/make_ofe_grid.md)
builds a regular grid (square cells, optional rotation and origin shift)
over the bounding box of the trial data and assigns each yield-monitor
point to a cell. A cell is **selected for analysis** if and only if:

- it contains at least `min_per_cell` observations, and
- all those observations come from a single treatment.

Mixed-treatment cells (those that straddle the boundary between strips)
are dropped. The per-cell response used downstream is the **median** of
the raw observations falling in that cell, which is robust to the
typical outliers of yield-monitor data.

*Implementation details:* the grid origin and cell size are passed
through deterministically; the same arguments always produce the same
cells, so `ofemt(grid = NULL, ...)` and
`ofemt(grid = make_ofe_grid(...))` are interchangeable.

## 3. Effective sample size (ESS)

A one-way ANOVA on cell medians yields residuals \\e_i\\. Their spatial
autocorrelation is estimated with the approximate profile likelihood
estimator (APLE) of \\\rho\\ (Li, Calder & Cressie, 2007), and the
**effective sample size** is computed from Griffith (2005):

\\ n^{\*} = n \left(1 - a \cdot \frac{n-1}{n} \cdot \left(1 -
\exp(-b\rho + c\sqrt{\rho})\right)\right) \\

with the coefficients given in the published formula. We also report
**Moran’s *I*** of the residuals as a sanity check on the magnitude of
the autocorrelation. \\\rho\\ is constrained to \\\[0, 1\]\\ before
entering the ESS formula.

## 4. Permutational ANOVA on subsamples

Given the ESS,
[`ofemt()`](https://ppaccioretti.github.io/ofemeantest/reference/ofemt.md)
runs `n_s` independent permutational ANOVAs. Each run:

1.  samples \\\lceil n^{\*}/k \rceil\\ cells **without replacement**
    from each of the \\k\\ treatments;
2.  fits a one-way ANOVA on the subsample;
3.  permutes the treatment labels `n_p` times and records the proportion
    of permuted *F*-statistics that exceed the observed one (this is the
    per-run *p*-value).

The reported test statistic for the comparison is the **median** of the
`n_s` per-run *p*-values. The full distribution is retained in
`res$perm_runs` and can be inspected with
[`plot_pvalue_hist()`](https://ppaccioretti.github.io/ofemeantest/reference/plot_pvalue_hist.md).

## 5. Multiple comparisons

When the trial contains more than two treatments, every pairwise
comparison is tested as above and the resulting median *p*-values are
adjusted for multiplicity via
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
(Bonferroni, Holm or Benjamini–Hochberg, selected by the
`p_adjust_method` argument). A compact letter display is then computed
from the adjusted *p*-values with
[`multcompView::multcompLetters()`](https://rdrr.io/pkg/multcompView/man/multcompLetters.html).

## 6. Assumptions and limitations

- The protocol assumes the trial layout produces enough cells per
  treatment for the subsampling step to be meaningful.
- The ESS formula is an approximation valid for approximately normally
  distributed data with non-negative autocorrelation.
- The current implementation supports a **single categorical
  treatment**. Multi-factor designs, factor interactions and continuous
  covariates are not yet supported.

## References

- Córdoba M., Paccioretti P., Balzarini M. (in review). A new method to
  compare treatments in unreplicated on-farm experimentation.
- Griffith, D.A. (2005). Effective geographic sample size in the
  presence of spatial autocorrelation. *Annals of the Association of
  American Geographers* 95(4): 740–760.
- Li, H., Calder, C.A. & Cressie, N. (2007). Beyond Moran’s *I*: Testing
  for spatial dependence based on the spatial autoregressive model.
  *Geographical Analysis* 39(4): 357–375.
- Vega A., Córdoba M., Balzarini M. (2019). Protocol for automating
  error removal from yield maps. *Precision Agriculture* 20: 1030–1044.
