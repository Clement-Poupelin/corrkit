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
