# Significance stars of p-values

`***` p \< 0.001, `**` p \< 0.01, `*` p \< 0.05; otherwise `ns` (or `.`
for p \< 0.10 when `trend = TRUE`). Missing p -\> "".

## Usage

``` r
p_to_stars(p, ns = "ns", trend = FALSE)
```

## Arguments

- p:

  numeric vector of p-values

- ns:

  text for a non-significant p

- trend:

  mark 0.05 \<= p \< 0.10 with "."

## Value

character vector, same length as `p`

## Examples

``` r
p_to_stars(c(0.0004, 0.003, 0.04, 0.07, 0.5, NA), ns = "", trend = TRUE)
#> [1] "***" "**"  "*"   "."   ""    ""   
```
