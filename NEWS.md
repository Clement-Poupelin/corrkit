# corrkit 0.2.0

* **One computing function, `cor_compute()`**, replaces `cor_test_matrix()` and `cor_cross()`:
  `y = NULL` gives every pair of `x` (matrix), `y` gives the cross analysis `x` x `y`, one `x`
  and one `y` a single pair; `covariates` works in every case (partial correlation matrices
  are new). One result class, `cor_result`.
* Methods of `cor_result`: `print()`, `summary()` (settings, variables, strongest pairs,
  content of the object), `plot()` and `autoplot()` with `type = "heatmap"`, `"forest"`
  or `"scatter"`, `as.matrix(what = )`, `as.data.frame()`.
* New `plot_cor_forest()`: each coefficient with its confidence interval.
* `plot_cor_pair()`: choice of the fitted curve (`smooth = "lm"`, `"loess"`, `"gam"`, `"glm"`,
  a function or `"none"`, plus `smooth_args` and `se`; default `"lm"` for Pearson and partial
  correlations, `"loess"` for Spearman and Kendall), colour-blind-safe default colours
  (Okabe-Ito: `point_colour`, `line_colour`, `group_colours`), shorter axis titles for the
  added-variable plot. The `palette` argument is removed.
* `tbl_cor(layout = "matrix" | "wide" | "long")` replaces `tbl_cor_matrix()`,
  `tbl_cor_cross()` and `tbl_cor_long()`.
* `plot_cor_heatmap()` takes a `cor_result` and draws a triangle (matrix) or a rectangle
  (cross); `plot_cor_cross()` is removed.
* `cor_binary_effects()` also accepts a `cor_result`.
* Removed: `cor_test_matrix()`, `cor_cross()`, `plot_cor_cross()`, `tbl_cor_matrix()`,
  `tbl_cor_cross()`, `tbl_cor_long()`, `partial_rank_cor()`.
* Plot texts are typeset: with the ggtext package (Suggests), subtitles, captions and axis titles
  are rendered as markdown, so that r_s, r_xy.z, e_x.z, r_rb... appear as an italic symbol with a
  subscript (`options(corrkit.markdown = FALSE)` for plain text); the colour legend uses plotmath.
* `options(corrkit.unicode = FALSE)` (default) writes the square root, minus sign and tau as
  `sqrt`, `-`, `tau`, so that figures can be saved with the pdf() device; `TRUE` keeps the symbols
  (PNG, or `ggsave(device = cairo_pdf)`).

# corrkit 0.1.0

* First version, built from `03_correlation.R`.
* Engines: `stats::cor.test()`, `ppcor::pcor.test()`, `psych::corr.test()`, `correlation::correlation()`
  (new engine, gives access to Kendall and to the robust methods of the correlation package).
* Confidence intervals (Fisher z, via `correlation::cor_to_ci()` when installed).
* New: `cor_cross()`, `plot_cor_pair()`, `cor_binary_effects()`, gt tables `tbl_cor_matrix()`,
  `tbl_cor_cross()`, `tbl_cor_long()`, `tbl_cor_binary()`, example data `simulate_cor_data()`.
* Three vignettes (in French).

## Changes after the first site

* p-values are now **raw by default** (`p_adjust = "none"`) in `cor_test_matrix()`, `cor_cross()`
  and `plot_cor_cross()`; `adjust_by` is used only when an adjustment is requested.
  `tbl_cor_long()` hides the adjusted column when p-values are raw.
* New vignettes: detailed workflow, multiple-testing corrections, comparison of correlation coefficients.
* Site now 100% Quarto, built with altdoc (`altdoc/`, workflow `.github/workflows/altdoc.yaml`);
  pkgdown removed. Vignettes in Quarto, with tabsets and callouts.
* Logo script in `data-raw/logo.R` (hexSticker).
