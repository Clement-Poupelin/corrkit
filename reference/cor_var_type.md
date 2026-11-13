# Type of a variable for a correlation

Type of a variable for a correlation

## Usage

``` r
cor_var_type(x)
```

## Arguments

- x:

  a vector (numeric, logical, factor or character)

## Value

`"binary"` (exactly 2 distinct non-missing values), `"ordinal"` (ordered
factor with \> 2 levels), `"nominal"` (unordered factor or character
with \> 2 levels) or `"continuous"`.

## Examples

``` r
cor_var_type(c(1.2, 3.4, 5.6))
#> [1] "continuous"
cor_var_type(factor(c("low", "high", "low")))
#> [1] "binary"
cor_var_type(factor(c("a", "b", "c")))
#> [1] "nominal"
```
