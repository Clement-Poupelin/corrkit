# Effect sizes of binary variables on continuous scores

For a binary variable (e.g. a biomarker "high" vs "low"), the
correlation is bounded by the split between groups and is not comparable
with the correlation of a continuous variable (see
[`cor_binary_note()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_note.md)).
This function returns what to report instead, with existing tests:

- difference in means M_1 - M_0 and its CI:
  [`stats::t.test()`](https://rdrr.io/r/stats/t.test.html) with equal
  variances (the test of the point-biserial r);

- rank-biserial correlation r_rb = 2U / (n_0 n_1) - 1 = P(Y_1 \> Y_0) -
  P(Y_1 \< Y_0) = 2 AUC - 1 (Kerby 2014), with U the Mann-Whitney
  statistic of
  [`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html); CI
  from
  [`effectsize::rank_biserial()`](https://easystats.github.io/effectsize/reference/rank_biserial.html)
  if installed;

- for comparison, the observed r_s and r_pb and the bound sqrt(3pq) of
  \|r_s\|.

## Usage

``` r
cor_binary_effects(data, bin_vars, y_vars, labels = NULL, ci = 0.95)
```

## Arguments

- data:

  data frame

- bin_vars:

  binary variables (non-binary ones are dropped with a warning)

- y_vars:

  continuous outcomes

- labels:

  named vector variable -\> display label (kept for the tables)

- ci:

  confidence level

## Value

tibble, one row per (binary variable, score): X, Y, level0, level1, n0,
n1, prop1, mean0, sd0, mean1, sd1, median0, median1, diff, diff_low,
diff_high, p_t, r_rb, r_rb_low, r_rb_high, p_mw, r_s, r_s_bound, r_pb.

## Details

Group 1 = upper level of the binary variable (as in
[`as_cor_numeric()`](https://clement-poupelin.github.io/corrkit/reference/as_cor_numeric.md)).
p-values are raw.

## References

Kerby DS (2014). Comprehensive Psychology 3:11.IT.3.1.
doi:10.2466/11.IT.3.1

## Examples

``` r
d <- simulate_cor_data(150)
cor_binary_effects(d, c("BIN_CRP", "BIN_IL6"), c("depression", "fatigue"))
#> # A tibble: 4 × 24
#>   X       Y      level0 level1    n0    n1 prop1 mean0   sd0 mean1   sd1 median0
#>   <chr>   <chr>  <chr>  <chr>  <int> <int> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#> 1 BIN_CRP depre… low    high      74    69 0.483 11.5   3.30 12.8   3.22    11.7
#> 2 BIN_CRP fatig… low    high      73    69 0.486  8.02  3.51  8.61  3.69     8.7
#> 3 BIN_IL6 depre… low    high     115    27 0.190 11.8   3.31 13.8   2.80    12  
#> 4 BIN_IL6 fatig… low    high     115    28 0.196  7.61  3.42 10.3   3.07     8  
#> # ℹ 12 more variables: median1 <dbl>, diff <dbl>, diff_low <dbl>,
#> #   diff_high <dbl>, p_t <dbl>, r_rb <dbl>, r_rb_low <dbl>, r_rb_high <dbl>,
#> #   p_mw <dbl>, r_s <dbl>, r_s_bound <dbl>, r_pb <dbl>
```
