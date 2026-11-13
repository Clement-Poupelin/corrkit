# Table of cross correlations X x Y (wide: one row per predictor, one column per score)

Same content and same colour scale as
[`plot_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_cross.md).
Binary predictors are grouped under "Binary 0/1" and carry a footnote.

## Usage

``` r
tbl_cor_cross(
  x,
  labels = NULL,
  digits = 2,
  palette = COR_PAL,
  show_note = TRUE,
  title = NULL,
  subtitle = NULL,
  stub_label = "Predictor"
)
```

## Arguments

- x:

  output of
  [`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md)
  or
  [`plot_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_cross.md)

- labels:

  named vector variable -\> label (default: those given to
  [`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md))

- digits:

  decimals of the coefficients

- palette:

  colours (named low / high); `NULL` = no colour

- show_note:

  add the reading lines of binary variables

- title, subtitle:

  override the automatic texts

- stub_label:

  header of the row labels

## Value

a `gt_tbl`

## Examples

``` r
d  <- simulate_cor_data(150)
cc <- cor_cross(d, c("CRP", "IL6", "BIN_CRP"), c("depression", "fatigue"),
                covariates = "age", labels = attr(d, "labels"))
tbl_cor_cross(cc)


  


Spearman correlations
```
