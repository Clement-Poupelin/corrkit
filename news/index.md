# Changelog

## corrkit 0.1.0

- First version, built from `03_correlation.R`.
- Engines: [`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html),
  [`ppcor::pcor.test()`](https://rdrr.io/pkg/ppcor/man/pcor.test.html),
  [`psych::corr.test()`](https://rdrr.io/pkg/psych/man/corr.test.html),
  [`correlation::correlation()`](https://easystats.github.io/correlation/reference/correlation.html)
  (new engine, gives access to Kendall and to the robust methods of the
  correlation package).
- Confidence intervals (Fisher z, via
  [`correlation::cor_to_ci()`](https://easystats.github.io/correlation/reference/cor_to_p.html)
  when installed).
- New:
  [`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md),
  [`plot_cor_pair()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_pair.md),
  [`cor_binary_effects()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_effects.md),
  gt tables
  [`tbl_cor_matrix()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_matrix.md),
  [`tbl_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_cross.md),
  [`tbl_cor_long()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_long.md),
  [`tbl_cor_binary()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_binary.md),
  example data
  [`simulate_cor_data()`](https://clement-poupelin.github.io/corrkit/reference/simulate_cor_data.md).
- Three vignettes (in French).
