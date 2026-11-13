# Table of a correlation matrix (lower triangle), same layout as the heatmap

Table of a correlation matrix (lower triangle), same layout as the
heatmap

## Usage

``` r
tbl_cor_matrix(
  cm,
  labels = NULL,
  use_adj = TRUE,
  order = c("original", "hclust"),
  digits = 2,
  palette = COR_PAL,
  bin_mark = " [0/1]",
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- cm:

  output of
  [`cor_test_matrix()`](https://clement-poupelin.github.io/corrkit/reference/cor_test_matrix.md)

- labels:

  named vector variable -\> label

- use_adj:

  stars on adjusted (`TRUE`) or raw p-values

- order:

  `"original"` or `"hclust"` (same ordering as
  [`plot_cor_heatmap()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_heatmap.md))

- digits:

  decimals

- palette:

  colours (named low / high); `NULL` = no colour

- bin_mark:

  marker of binary variables

- title, subtitle:

  override the automatic texts

## Value

a `gt_tbl`

## Examples

``` r
d  <- simulate_cor_data(150)
cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"))
tbl_cor_matrix(cm, labels = attr(d, "labels"))


  


Correlation matrix
```
