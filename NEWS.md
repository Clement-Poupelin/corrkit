# corrkit 0.1.0

* First version, built from `03_correlation.R`.
* Engines: `stats::cor.test()`, `ppcor::pcor.test()`, `psych::corr.test()`, `correlation::correlation()`
  (new engine, gives access to Kendall and to the robust methods of the correlation package).
* Confidence intervals (Fisher z, via `correlation::cor_to_ci()` when installed).
* New: `cor_cross()`, `plot_cor_pair()`, `cor_binary_effects()`, gt tables `tbl_cor_matrix()`,
  `tbl_cor_cross()`, `tbl_cor_long()`, `tbl_cor_binary()`, example data `simulate_cor_data()`.
* Three vignettes (in French).
