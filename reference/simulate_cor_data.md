# Simulated biomarkers, clinical scores and covariates

Data generated on the fly (no stored data set), with known structure:

- a latent inflammation L drives `CRP` (log-normal, skewed), `IL6`,
  `BDNF` (negatively) and the `depression` score;

- `IL6` drives `fatigue`;

- `age` drives both `NfL` and `cognition`: their raw correlation is
  mostly confounded by age and shrinks once adjusted for age;

- `BIN_CRP` = CRP above its median ("high" vs "low", 50/50 split) and
  `BIN_IL6` = IL6 above its 80th percentile (20/80 split): the same
  information as the continuous markers, cut into two groups;

- about `na_rate` missing values in the biomarkers and scores.

## Usage

``` r
simulate_cor_data(n = 200, seed = 2026, na_rate = 0.03)
```

## Arguments

- n:

  number of individuals

- seed:

  seed of the simulation

- na_rate:

  share of missing values in each biomarker and score

## Value

data.frame with an attribute `"labels"` (named vector variable -\>
label)

## Details

The random-number state of the session is restored on exit.

## Examples

``` r
d <- simulate_cor_data(200)
str(d)
#> 'data.frame':    200 obs. of  11 variables:
#>  $ age       : num  61 42 57 54 47 25 46 43 56 49 ...
#>  $ sex       : Factor w/ 2 levels "F","M": 1 1 1 2 2 2 1 1 2 2 ...
#>  $ CRP       : num  1.95 0.98 0.51 3.89 0.55 0.43 0.6 1.66 1.46 2.39 ...
#>  $ IL6       : num  0.82 -0.65 -0.28 0.71 -0.07 -0.45 -0.29 1.5 -1.36 0.14 ...
#>  $ BDNF      : num  0.3 -0.35 -0.21 0.89 1.55 1.31 0.32 0.17 -1.48 -1.3 ...
#>  $ NfL       : num  0.93 -0.15 -0.76 -0.86 -1.08 -1.47 -1.5 NA -0.99 -0.06 ...
#>  $ depression: num  7.3 7.3 NA 14.4 9.1 10.2 6 14.1 15.2 11.4 ...
#>  $ fatigue   : num  4 5.6 4.6 9.2 9.1 7.6 7.5 12.3 4.5 0.1 ...
#>  $ cognition : num  23.6 25.7 25.5 24.3 26 NA 27.6 26.3 26.1 27.7 ...
#>  $ BIN_CRP   : Factor w/ 2 levels "low","high": 2 1 1 2 1 1 1 2 2 2 ...
#>  $ BIN_IL6   : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 2 1 1 ...
#>  - attr(*, "labels")= Named chr [1:11] "Age" "Sex" "CRP" "IL-6" ...
#>   ..- attr(*, "names")= chr [1:11] "age" "sex" "CRP" "IL6" ...
attr(d, "labels")
#>                age                sex                CRP                IL6 
#>              "Age"              "Sex"              "CRP"             "IL-6" 
#>               BDNF                NfL         depression            fatigue 
#>             "BDNF"              "NfL" "Depression score"    "Fatigue score" 
#>          cognition            BIN_CRP            BIN_IL6 
#>  "Cognition score"         "CRP high"        "IL-6 high" 
```
