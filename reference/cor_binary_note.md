# Reading note for a correlation involving a binary (0/1) variable

- **Pearson** + binary = point-biserial r_pb. Its test is Student's
  t-test (equal variances); partial (adjusted) = slope test of
  `lm(Y ~ X + covariates)` (Frisch-Waugh-Lovell: same t, same df).
  \|r_pb\| depends on the split p/q: for a biomarker categorised from a
  continuous measurement, r_pb = r x phi(c) / sqrt(p q) (c = cut-off in
  SD units), i.e. at most 0.80 x r at the median, 0.70 x r for a 20/80
  split.

- **Spearman** + binary = Pearson on ranks; a linear function of the
  Mann-Whitney U (same test, approximately). r_s = r_rb x sqrt(3 p q)
  (no ties in Y), with r_rb = 2 AUC - 1 the rank-biserial correlation:
  \|r_s\| \<= 0.87 (50/50), 0.69 (20/80).

- **Kendall** + binary: the numerator of tau (concordant - discordant
  pairs) is n_0 n_1 r_rb, i.e. the Mann-Whitney comparison; its maximum
  depends on the split and on the ties.

## Usage

``` r
cor_binary_note(
  meth = c("spearman", "pearson", "kendall"),
  lang = c("en", "fr"),
  short = FALSE,
  mark = " [0/1]"
)
```

## Arguments

- meth:

  `"spearman"`, `"pearson"` or `"kendall"`

- lang:

  `"en"` (figure captions, table notes) or `"fr"` (callouts of a French
  report)

- short:

  `TRUE` = one line for a caption; `FALSE` = full note with formulas and
  references

- mark:

  marker of binary variables used in the text (same as `bin_mark` of the
  plots)

## Value

a character string

## Details

In every case the effect size to report is r_rb or the difference
between groups
([`cor_binary_effects()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_effects.md)),
not the correlation.

## References

Kerby DS (2014). The simple difference formula: an approach to teaching
nonparametric correlation. *Comprehensive Psychology* 3:11.IT.3.1.
doi:10.2466/11.IT.3.1

Babchishin KM, Helmus LM (2016). The influence of base rates on
correlations. *Behavior Research Methods* 48:1021-1031.
doi:10.3758/s13428-015-0627-7

## Examples

``` r
cor_binary_note("spearman", short = TRUE)
#> [1] "[0/1] binary: r_s ≈ Mann-Whitney (rank difference between groups); |r_s| ≤ √(3pq) ≤ 0.87, not comparable with continuous variables"
cat(cor_binary_note("pearson", lang = "fr"))
#> Variables binaires (notées [0/1], classe haute = 1) : le r de Pearson est la corrélation point-bisériale. Son test est le test t de Student entre les deux groupes (variances égales) et, ajusté sur des covariables, le test de la pente de lm(score ~ groupe + covariables) : même t, mêmes degrés de liberté. Sa valeur dépend de la répartition p/q entre groupes : pour un biomarqueur catégorisé à partir d'une mesure continue, r_pb = r × φ(c) / √(pq) (c = seuil en unités d'écart-type), soit au plus 0,80 × r à la médiane et 0,70 × r pour une répartition 20/80. Il n'est donc pas comparable au r d'un biomarqueur continu : la taille de l'effet se lit sur la différence de score moyen entre groupes. (Kerby 2014, Comprehensive Psychology, doi:10.2466/11.IT.3.1; Babchishin & Helmus 2016, Behav Res Methods, doi:10.3758/s13428-015-0627-7).
```
