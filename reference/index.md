# Package index

## Calculer

Chaque fonction appelle un moteur existant (stats, ppcor, psych,
correlation) et ajuste les p-values sur une famille explicite.

- [`cor_test_matrix()`](https://clement-poupelin.github.io/corrkit/reference/cor_test_matrix.md)
  : Correlation and p-value matrices of all pairs (pairwise complete
  observations)
- [`cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/cor_cross.md)
  : Cross correlations X x Y (e.g. biomarkers x scores), optionally
  partial
- [`partial_cor_test()`](https://clement-poupelin.github.io/corrkit/reference/partial_cor_test.md)
  [`partial_rank_cor()`](https://clement-poupelin.github.io/corrkit/reference/partial_cor_test.md)
  : Correlation test between x and y, partial if covariates are given
- [`cor_binary_effects()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_effects.md)
  : Effect sizes of binary variables on continuous scores

## Graphiques

Objets ggplot, avec sous-titre (méthode, moteur, N) et légende de
lecture.

- [`plot_cor_heatmap()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_heatmap.md)
  : Lower-triangular heatmap of a correlation matrix
- [`plot_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_cross.md)
  : Heatmap of cross correlations X x Y, optionally partial
- [`plot_cor_pair()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_pair.md)
  : Plot of one correlation, with its coefficient, CI and test in the
  subtitle

## Tableaux

Tableaux gt, même échelle de couleurs que les heatmaps, notes de bas de
tableau.

- [`tbl_cor_binary()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_binary.md)
  : Table of the effect sizes of binary variables
- [`tbl_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_cross.md)
  : Table of cross correlations X x Y (wide: one row per predictor, one
  column per score)
- [`tbl_cor_long()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_long.md)
  : Long table of correlations: one row per pair, with CI, raw and
  adjusted p
- [`tbl_cor_matrix()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_matrix.md)
  : Table of a correlation matrix (lower triangle), same layout as the
  heatmap

## Notes de lecture et utilitaires

- [`cor_binary_note()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_note.md)
  : Reading note for a correlation involving a binary (0/1) variable

- [`cor_var_type()`](https://clement-poupelin.github.io/corrkit/reference/cor_var_type.md)
  : Type of a variable for a correlation

- [`as_cor_numeric()`](https://clement-poupelin.github.io/corrkit/reference/as_cor_numeric.md)
  : Numeric coding of a variable for a correlation

- [`mark_binary()`](https://clement-poupelin.github.io/corrkit/reference/mark_binary.md)
  :

  Add a marker (default `" [0/1]"`) to the labels of binary variables

- [`p_to_stars()`](https://clement-poupelin.github.io/corrkit/reference/p_to_stars.md)
  : Significance stars of p-values

- [`theme_pub()`](https://clement-poupelin.github.io/corrkit/reference/theme_pub.md)
  : ggplot2 theme of the corrkit figures

- [`COR_PAL`](https://clement-poupelin.github.io/corrkit/reference/COR_PAL.md)
  : Palette of the correlation heatmaps and tables (named low / high,
  white at 0)

## Données d’exemple

- [`simulate_cor_data()`](https://clement-poupelin.github.io/corrkit/reference/simulate_cor_data.md)
  : Simulated biomarkers, clinical scores and covariates
