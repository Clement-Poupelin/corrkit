# Long table of correlations: one row per pair, with CI, raw and adjusted p

Works for
[`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md)
/
[`plot_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_cross.md)
outputs (rows grouped by the family of the adjustment) and for
[`cor_test_matrix()`](https://clement-poupelin.github.io/corrkit/reference/cor_test_matrix.md)
outputs (all pairs). Rows with an adjusted p below `alpha` are in bold.

## Usage

``` r
tbl_cor_long(
  x,
  labels = NULL,
  digits = 2,
  alpha = 0.05,
  palette = COR_PAL,
  title = NULL,
  subtitle = NULL,
  sort_by = c("none", "abs")
)
```

## Arguments

- x:

  output of
  [`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md),
  [`plot_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_cross.md)
  or
  [`cor_test_matrix()`](https://clement-poupelin.github.io/corrkit/reference/cor_test_matrix.md)

- labels:

  named vector variable -\> label

- digits:

  decimals of the coefficients

- alpha:

  threshold for bold rows

- palette:

  colours (named low / high) of the coefficient column; `NULL` = no
  colour

- title, subtitle:

  override the automatic texts

- sort_by:

  `"none"` (order of the variables) or `"abs"` (decreasing
  \|coefficient\| within groups)

## Value

a `gt_tbl`

## Examples

``` r
d  <- simulate_cor_data(150)
cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "BIN_CRP"))
tbl_cor_long(cm, labels = attr(d, "labels"))


  


Spearman correlations
```
