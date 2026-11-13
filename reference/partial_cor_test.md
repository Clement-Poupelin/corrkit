# Correlation test between x and y, partial if covariates are given

- no covariate:
  [`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html)
  (`exact = FALSE`);

- covariates:
  [`ppcor::pcor.test()`](https://rdrr.io/pkg/ppcor/man/pcor.test.html)
  (default engine); fallback or `engine = "base"`: correlation of the
  residuals of x and y on Z (for Spearman x, y AND Z are ranked, as in
  ppcor), t-test with n - 2 - k df. For Pearson this is exactly the
  slope test of `lm(y ~ x + Z)` (Frisch-Waugh-Lovell theorem). Kendall
  partial correlations need ppcor.

`partial_rank_cor()` is the former name, kept for backward
compatibility.

## Usage

``` r
partial_cor_test(
  x,
  y,
  Z = NULL,
  meth = c("spearman", "pearson", "kendall"),
  engine = c("ppcor", "base"),
  ci = 0.95
)

partial_rank_cor(x, y, Z = NULL, meth = c("spearman", "pearson"))
```

## Arguments

- x, y:

  numeric vectors (binary: 0/1, see
  [`as_cor_numeric()`](https://clement-poupelin.github.io/corrkit/reference/as_cor_numeric.md))

- Z:

  data frame of covariates (factors allowed: dummy coded) or `NULL`

- meth:

  `"spearman"` (ranks), `"pearson"` (values) or `"kendall"`

- engine:

  `"ppcor"` (default) or `"base"` (residual method)

- ci:

  confidence level of the interval (`NULL` = none); Fisher z on n - k

## Value

list(estimate, p.value, n, k, ci_low, ci_high, engine, ci_engine)

## Examples

``` r
d <- simulate_cor_data(100)
partial_cor_test(d$IL6, d$fatigue, meth = "spearman")
#> $estimate
#> [1] 0.5246237
#> 
#> $p.value
#> [1] 5.712139e-08
#> 
#> $n
#> [1] 94
#> 
#> $k
#> [1] 0
#> 
#> $ci_low
#> [1] 0.3550094
#> 
#> $ci_high
#> [1] 0.6607999
#> 
#> $engine
#> [1] "stats::cor.test"
#> 
#> $ci_engine
#> [1] "correlation::cor_to_ci"
#> 
partial_cor_test(d$NfL, d$cognition, Z = d[, c("age", "sex")], meth = "pearson")
#> $estimate
#> [1] -0.04689676
#> 
#> $p.value
#> [1] 0.6571061
#> 
#> $n
#> [1] 94
#> 
#> $k
#> [1] 2
#> 
#> $ci_low
#> [1] -0.2493194
#> 
#> $ci_high
#> [1] 0.1594522
#> 
#> $engine
#> [1] "ppcor::pcor.test"
#> 
#> $ci_engine
#> [1] "correlation::cor_to_ci"
#> 
```
