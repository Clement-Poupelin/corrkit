# Lower-triangular heatmap of a correlation matrix

Lower-triangular heatmap of a correlation matrix

## Usage

``` r
plot_cor_heatmap(
  cm,
  labels = NULL,
  use_adj = TRUE,
  order = c("hclust", "original"),
  show_values = TRUE,
  base_size = 10,
  title = NULL,
  value_size = 2.6,
  show_note = TRUE,
  bin_mark = " [0/1]",
  split_binary = TRUE,
  palette = COR_PAL,
  caption_width = 110,
  separator = c("facet", "line")
)
```

## Arguments

- cm:

  output of
  [`cor_test_matrix()`](https://clement-poupelin.github.io/corrkit/reference/cor_test_matrix.md)

- labels:

  named vector variable -\> display label

- use_adj:

  stars on adjusted (`TRUE`) or raw p-values

- order:

  `"hclust"` (within each block if `split_binary`) or `"original"`

- show_values:

  print the coefficient and the stars in each cell

- base_size:

  base font size

- title:

  title (default "Correlation matrix")

- value_size:

  size of the values

- show_note:

  add the binary-variable line to the caption

- bin_mark:

  marker appended to the labels of binary variables

- split_binary:

  binary variables in a separate block (last rows / columns)

- palette:

  colours of the negative / positive correlations (named low / high)

- caption_width:

  wrap the caption lines at this number of characters (`NULL` = no
  wrapping)

- separator:

  `"facet"` (2 x 2 panels with strips, default) or `"line"` (separator
  lines, square tiles)

## Value

a ggplot

## Examples

``` r
d  <- simulate_cor_data(150)
cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "NfL", "depression", "BIN_CRP", "BIN_IL6"))
plot_cor_heatmap(cm, labels = attr(d, "labels"))
```
