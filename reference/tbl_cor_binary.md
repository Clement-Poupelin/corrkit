# Table of the effect sizes of binary variables

Formats the output of
[`cor_binary_effects()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_effects.md):
groups and sizes, means (SD), difference in means with its CI and
t-test, rank-biserial correlation r_rb with its CI and Mann-Whitney
test, and the observed r_s next to its maximum for this split, so that
the attenuation of the correlation is visible.

## Usage

``` r
tbl_cor_binary(eff, labels = NULL, digits = 2, title = NULL, subtitle = NULL)
```

## Arguments

- eff:

  output of
  [`cor_binary_effects()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_effects.md)

- labels:

  named vector variable -\> label

- digits:

  decimals

- title, subtitle:

  override the automatic texts

## Value

a `gt_tbl`

## Examples

``` r
d <- simulate_cor_data(150)
tbl_cor_binary(cor_binary_effects(d, c("BIN_CRP", "BIN_IL6"), c("depression", "fatigue")),
               labels = attr(d, "labels"))


  


Binary variables: effect sizes to report
```
