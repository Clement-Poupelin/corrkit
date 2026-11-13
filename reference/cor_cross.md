# Cross correlations X x Y (e.g. biomarkers x scores), optionally partial

Every pair (x, y) is tested with
[`partial_cor_test()`](https://clement-poupelin.github.io/corrkit/reference/partial_cor_test.md);
p-values are then adjusted with
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) within the
family chosen by `adjust_by`. The result feeds
[`plot_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_cross.md),
[`tbl_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_cross.md)
and
[`tbl_cor_long()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_long.md).

## Usage

``` r
cor_cross(
  data,
  x_vars,
  y_vars,
  covariates = NULL,
  labels = NULL,
  p_adjust = "BH",
  adjust_by = c("y", "x", "all", "none"),
  meth = c("spearman", "pearson", "kendall"),
  engine = c("ppcor", "base"),
  ci = 0.95,
  bin_mark = " [0/1]"
)
```

## Arguments

- data:

  data frame

- x_vars:

  predictors (rows of the heatmap, e.g. biomarkers)

- y_vars:

  outcomes (columns, e.g. scores)

- covariates:

  covariates of the partial correlation (`NULL` = unadjusted)

- labels:

  named vector variable -\> display label

- p_adjust:

  method of [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)

- adjust_by:

  family of the adjustment: `"y"` (the predictors of each score,
  default, as in univariable forest plots), `"x"` (the scores of each
  predictor), `"all"` (every test), `"none"`

- meth:

  `"spearman"` (default), `"pearson"` or `"kendall"`

- engine:

  `"ppcor"` (default) or `"base"` for the partial correlations

- ci:

  confidence level (`NULL` = none)

- bin_mark:

  marker appended to the labels of binary variables

## Value

object of class `cor_cross`: list(table, meta). `table` has one row per
pair: X, Y, x_type, y_type, rho (the coefficient, whatever the method),
ci_low, ci_high, p, n, k, engine, p_adj.

## Examples

``` r
d  <- simulate_cor_data(150)
cc <- cor_cross(d, x_vars = c("CRP", "IL6", "NfL", "BIN_CRP"),
                y_vars = c("depression", "fatigue", "cognition"),
                covariates = c("age", "sex"), labels = attr(d, "labels"))
cc
#> <cor_cross> 4 x 3 Spearman correlations (partial, adjusted for age, sex)
#>   Adjusted p, Benjamini-Hochberg (FDR) over the 4 predictors of each score (column): * < 0.05, ** < 0.01, *** < 0.001; no star = adjusted p ≥ 0.05.
#>        X depression fatigue cognition
#>      CRP      0.22*    0.07     -0.01
#>      IL6    0.34*** 0.38***      0.00
#>      NfL       0.02    0.18    -0.25*
#>  BIN_CRP      0.18*    0.08      0.02
```
