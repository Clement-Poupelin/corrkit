# corrkit

Analyses de corrélation avec des heatmaps et des tableaux qui expliquent ce qu'ils montrent.

`corrkit` ne recalcule rien : chaque coefficient, chaque p-value et chaque intervalle vient d'une fonction existante et validée (`stats::cor.test()`, `ppcor::pcor.test()`, `psych::corr.test()`, `correlation::correlation()`, `stats::p.adjust()`, `effectsize::rank_biserial()`). Le package ajoute par-dessus :

* **une famille de tests explicite** pour l'ajustement des p-values (`adjust_by`, Benjamini-Hochberg par défaut) ;
* **un traitement à part des variables binaires** (biomarqueur haut/bas, sexe) : blocs séparés, marque `[0/1]`, note de lecture, tailles d'effet à rapporter ;
* **des sorties annotées** : sous-titres, légendes et notes de tableau qui donnent la méthode, le moteur, les formules (notation de Makowski et al. 2020, [doi:10.21105/joss.02306](https://doi.org/10.21105/joss.02306)), l'ajustement et la lecture.

## Installation

```r
# depuis le dossier du package
install.packages(c("dplyr", "ggplot2", "gt", "purrr", "tibble", "tidyr"))
install.packages(c("ppcor", "psych", "correlation", "effectsize"))   # recommandés
devtools::document()      # génère man/ à partir des commentaires roxygen
devtools::install(build_vignettes = TRUE)
```

## En bref

```r
library(corrkit)
d   <- simulate_cor_data(200)
lab <- attr(d, "labels")

# Toutes les paires
cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"), method = "spearman")
plot_cor_heatmap(cm, labels = lab)
tbl_cor_matrix(cm, labels = lab)
tbl_cor_long(cm, labels = lab)

# Biomarqueurs x scores, ajustés sur l'âge et le sexe
cc <- cor_cross(d, c("CRP", "IL6", "NfL", "BIN_CRP"), c("depression", "fatigue", "cognition"),
                covariates = c("age", "sex"), labels = lab)
plot_cor_cross(cc)$plot
tbl_cor_cross(cc)
tbl_cor_long(cc)

# Une paire, et ce que mesure son coefficient
plot_cor_pair(d, "NfL", "cognition", covariates = "age", meth = "pearson", labels = lab)

# Variables binaires : ce qu'il faut rapporter
tbl_cor_binary(cor_binary_effects(d, "BIN_CRP", c("depression", "fatigue")), labels = lab)
```

## Fonctions

| Famille | Fonctions |
|---|---|
| Calcul | `cor_test_matrix()`, `cor_cross()`, `partial_cor_test()`, `cor_binary_effects()` |
| Graphiques (`ggplot`) | `plot_cor_heatmap()`, `plot_cor_cross()`, `plot_cor_pair()` |
| Tableaux (`gt`) | `tbl_cor_matrix()`, `tbl_cor_cross()`, `tbl_cor_long()`, `tbl_cor_binary()` |
| Notes et outils | `cor_binary_note()`, `cor_var_type()`, `as_cor_numeric()`, `mark_binary()`, `p_to_stars()`, `theme_pub()`, `COR_PAL` |
| Données | `simulate_cor_data()` |

## Vignettes

* `vignette("corrkit")` : prise en main, moteurs, tests multiples ;
* `vignette("partial-cross")` : corrélations croisées et partielles, familles d'ajustement, équivalences ppcor / résidus / `lm()` ;
* `vignette("binary-variables")` : variables binaires, atténuation, r rang-bisérial.

## Depuis `03_correlation.R`

Toutes les fonctions du fichier d'origine existent, avec des signatures compatibles (nouveaux arguments ajoutés en fin de liste) (`cor_var_type()`, `as_cor_numeric()`, `mark_binary()`, `cor_binary_note()`, `partial_cor_test()`, `partial_rank_cor()`, `cor_test_matrix()`, `plot_cor_heatmap()`, `plot_cor_cross()`, `COR_PAL`). Ajouts : argument `ci` et colonnes `ci_low` / `ci_high`, méthode Kendall, moteur `correlation`, classes `cor_matrix` / `cor_cross` avec `print()`, et `plot_cor_cross()` accepte aussi un objet `cor_cross`. Les helpers externes (`keep_existing()`, `get_label()`, `p_to_stars()`, `theme_pub()`, `rotate_x()`, `%||%`) sont maintenant internes au package.
