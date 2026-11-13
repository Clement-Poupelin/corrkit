# Numeric coding of a variable for a correlation

- logical -\> 0/1;

- binary factor or character -\> 0/1, the upper level coded 1 (so r \> 0
  means "higher in the upper level"); binary numeric -\> 0/1 (largest
  value = 1);

- ordered factor with \> 2 levels -\> 0, 1, 2, ... (Pearson then assumes
  equal spacing between levels: warning);

- unordered factor with \> 2 levels -\> NA with a warning (no meaningful
  correlation).

## Usage

``` r
as_cor_numeric(x, name = "x", meth = "spearman")
```

## Arguments

- x:

  vector to code

- name:

  name used in the warnings

- meth:

  correlation method (`"pearson"` triggers the equal-spacing warning)

## Value

numeric vector

## Details

Recoding a binary variable to 0/1 is a positive linear transformation:
it changes neither r, nor r_s, nor tau, nor their p-values.

## Examples

``` r
as_cor_numeric(factor(c("low", "high", "low"), levels = c("low", "high")))
#> [1] 0 1 0
```
