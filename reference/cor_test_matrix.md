# Correlation and p-value matrices of all pairs (pairwise complete observations)

The engine always returns **raw** p-values; corrkit then adjusts them
with [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) over
the k(k - 1)/2 distinct pairs, as
[`psych::corr.test()`](https://rdrr.io/pkg/psych/man/corr.test.html)
does above its diagonal.

## Usage

``` r
cor_test_matrix(
  data,
  vars,
  method = "spearman",
  p_adjust = "BH",
  engine = c("auto", "psych", "correlation", "stats"),
  ci = 0.95
)
```

## Arguments

- data:

  data frame

- vars:

  variables to correlate

- method:

  correlation method (see above)

- p_adjust:

  method of [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
  (`"BH"` by default, `"none"` for raw p)

- engine:

  `"auto"`, `"psych"`, `"correlation"` or `"stats"`

- ci:

  confidence level of the intervals (`NULL` = none)

## Value

object of class `cor_matrix`: list(r, p, p_adj, n, ci_low, ci_high,
long, types, method, p_adjust, engine, ci, ci_engine); `long` = one row
per pair (i \< j) with Var1, Var2, type1, type2, r, ci_low, ci_high, p,
p_adj, n.

## Details

Engines (`engine = "auto"` picks the first available in this order):

- Pearson / Spearman: `psych::corr.test(adjust = "none")`, then
  `correlation::correlation(p_adjust = "none")`, then
  [`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html) pair by
  pair;

- Kendall:
  [`correlation::correlation()`](https://easystats.github.io/correlation/reference/correlation.html),
  then [`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html)
  (psych tests tau with the t formula of r, so it is not used for
  Kendall);

- other methods of the correlation package (`"biweight"`,
  `"percentage"`, `"blomqvist"`, `"hoeffding"`, `"gaussian"`,
  `"shepherd"`, `"distance"`):
  [`correlation::correlation()`](https://easystats.github.io/correlation/reference/correlation.html)
  only.

Binary variables are coded 0/1 (see
[`as_cor_numeric()`](https://clement-poupelin.github.io/corrkit/reference/as_cor_numeric.md))
and flagged in `types`.

## Examples

``` r
d  <- simulate_cor_data(120)
cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"), method = "spearman")
cm
#> <cor_matrix> Spearman correlations, 5 variables (1 binary), 10 pairs
#>   engine: psych::corr.test() {psych} | p adjustment: Benjamini-Hochberg (FDR) over the 10 pairs | N = 112–116
#>         CRP     IL6    BDNF  NfL  BIN_CRP
#> CRP     1                                
#> IL6     0.36*** 1                        
#> BDNF    -0.16   -0.10  1                 
#> NfL     0.18    0.01   0.19  1           
#> BIN_CRP 0.87*** 0.33** -0.11 0.17 1      
head(cm$long)
#> # A tibble: 6 × 10
#>   Var1  Var2  type1      type2          r   ci_low ci_high       p     n   p_adj
#>   <chr> <chr> <chr>      <chr>      <dbl>    <dbl>   <dbl>   <dbl> <dbl>   <dbl>
#> 1 CRP   IL6   continuous contin…  0.356    0.177    0.512  1.18e-4   112 5.90e-4
#> 2 CRP   BDNF  continuous contin… -0.162   -0.342    0.0302 8.87e-2   112 1.27e-1
#> 3 IL6   BDNF  continuous contin… -0.0955  -0.281    0.0972 3.17e-1   112 3.52e-1
#> 4 CRP   NfL   continuous contin…  0.177   -0.0149   0.355  6.27e-2   112 1.25e-1
#> 5 IL6   NfL   continuous contin…  0.00784 -0.183    0.198  9.35e-1   112 9.35e-1
#> 6 BDNF  NfL   continuous contin…  0.194    0.00327  0.371  4.03e-2   112 1.01e-1
```
