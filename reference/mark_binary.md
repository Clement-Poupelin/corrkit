# Add a marker (default `" [0/1]"`) to the labels of binary variables

Add a marker (default `" [0/1]"`) to the labels of binary variables

## Usage

``` r
mark_binary(lbl, types, mark = " [0/1]")
```

## Arguments

- lbl:

  labels

- types:

  types from
  [`cor_var_type()`](https://clement-poupelin.github.io/corrkit/reference/cor_var_type.md),
  same length as `lbl`

- mark:

  marker appended to the binary labels (explained in captions and notes)
