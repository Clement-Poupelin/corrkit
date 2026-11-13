# Heatmap of cross correlations X x Y, optionally partial

Either computes the correlations
([`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md))
or draws an existing `cor_cross` object passed as `data`.

## Usage

``` r
plot_cor_cross(
  data,
  x_vars,
  y_vars,
  covariates = NULL,
  labels = NULL,
  p_adjust = "BH",
  adjust_by = c("y", "x", "all", "none"),
  meth = c("spearman", "pearson", "kendall"),
  engine = c("ppcor", "base"),
  base_size = 10,
  title = NULL,
  value_size = 2.8,
  show_note = TRUE,
  bin_mark = " [0/1]",
  split_binary = TRUE,
  palette = COR_PAL,
  caption_width = 90,
  ci = 0.95
)
```

## Arguments

- data:

  data frame, or the output of
  [`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md)

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

- base_size:

  base font size

- title:

  title (default " correlations")

- value_size:

  size of the values

- show_note:

  add the short reading note for binary variables to the caption

- bin_mark:

  marker appended to the labels of binary variables

- split_binary:

  binary variables in a separate facet ("Continuous" / "Binary 0/1")

- palette:

  colours of the negative / positive correlations (named low / high)

- caption_width:

  wrap the caption lines at this number of characters (`NULL` = no
  wrapping); ~60 for two heatmaps side by side

- ci:

  confidence level (`NULL` = none)

## Value

list(table, plot, note): `table` = the table of
[`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md)
(usable by
[`tbl_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_cross.md)
and
[`tbl_cor_long()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_long.md));
`note` = full reading note for binary variables (`NULL` if none)

## Examples

``` r
d <- simulate_cor_data(150)
res <- plot_cor_cross(d, x_vars = c("CRP", "IL6", "NfL", "BIN_CRP", "BIN_IL6"),
                      y_vars = c("depression", "fatigue", "cognition"),
                      covariates = c("age", "sex"), labels = attr(d, "labels"))
res$plot
```
