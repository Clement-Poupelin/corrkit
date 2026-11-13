# corrkit

**corrkit** réalise des analyses de corrélation avec des heatmaps et des
tableaux qui expliquent ce qu’ils montrent.

Le package ne recalcule rien : chaque coefficient, chaque p-value et
chaque intervalle vient d’une fonction existante et validée
([`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html),
[`ppcor::pcor.test()`](https://rdrr.io/pkg/ppcor/man/pcor.test.html),
[`psych::corr.test()`](https://rdrr.io/pkg/psych/man/corr.test.html),
[`correlation::correlation()`](https://easystats.github.io/correlation/reference/correlation.html),
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html),
[`effectsize::rank_biserial()`](https://easystats.github.io/effectsize/reference/rank_biserial.html)).
Il ajoute par-dessus :

- **une famille de tests explicite** pour l’ajustement des p-values
  (`adjust_by`, Benjamini-Hochberg par défaut) ;
- **un traitement à part des variables binaires** (biomarqueur haut/bas,
  sexe) : blocs séparés, marque `[0/1]`, note de lecture, tailles
  d’effet à rapporter ;
- **des sorties annotées** : sous-titres, légendes et notes de tableau
  qui donnent la méthode, le moteur, les formules (notation de Makowski
  et al. 2020,
  [doi:10.21105/joss.02306](https://doi.org/10.21105/joss.02306)),
  l’ajustement et la lecture.

## Installation

Version de développement, depuis GitHub :

``` r

# install.packages("pak")
pak::pak("clement-poupelin/corrkit")
```

Packages recommandés (moteurs et intervalles ; `corrkit` a un repli en R
de base sans eux) :

``` r

install.packages(c("ppcor", "psych", "correlation", "effectsize"))
```

## Données d’exemple

[`simulate_cor_data()`](https://clement-poupelin.github.io/corrkit/reference/simulate_cor_data.md)
génère des biomarqueurs, des scores cliniques et des covariables dont on
connaît la structure (effet confondant de l’âge, versions binaires des
marqueurs coupées à 50/50 et 20/80).

``` r

library(corrkit)
d   <- simulate_cor_data(200)
lab <- attr(d, "labels")      # libellés des variables
```

## Démarrage rapide

**Toutes les paires** d’un ensemble de variables :

``` r

cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"), method = "spearman")
plot_cor_heatmap(cm, labels = lab)   # heatmap annotée (ggplot)
tbl_cor_matrix(cm, labels = lab)     # même résultat en tableau (gt)
tbl_cor_long(cm, labels = lab)       # une ligne par paire : IC, p brute, p ajustée
```

**Biomarqueurs × scores**, ajustés sur l’âge et le sexe :

``` r

cc <- cor_cross(d, c("CRP", "IL6", "NfL", "BIN_CRP"), c("depression", "fatigue", "cognition"),
                covariates = c("age", "sex"), labels = lab)
plot_cor_cross(cc)$plot
tbl_cor_cross(cc)
```

**Une paire**, et ce que mesure son coefficient :

``` r

plot_cor_pair(d, "NfL", "cognition", covariates = "age", meth = "pearson", labels = lab)
```

**Variables binaires** : ce qu’il faut rapporter à la place de la
corrélation :

``` r

tbl_cor_binary(cor_binary_effects(d, "BIN_CRP", c("depression", "fatigue")), labels = lab)
```

## Pour aller plus loin

- [Prise en
  main](https://clement-poupelin.github.io/corrkit/articles/corrkit.html)
  : moteurs, heatmaps, tableaux, tests multiples ;
- [Corrélations croisées et
  partielles](https://clement-poupelin.github.io/corrkit/articles/partial-cross.html)
  : familles d’ajustement, équivalences ppcor / résidus /
  [`lm()`](https://rdrr.io/r/stats/lm.html) ;
- [Variables
  binaires](https://clement-poupelin.github.io/corrkit/articles/binary-variables.html)
  : atténuation, r rang-bisérial, ce qu’il faut rapporter ;
- [Référence des
  fonctions](https://clement-poupelin.github.io/corrkit/reference/index.html).

## Citation

``` r

citation("corrkit")
```

Citez aussi les moteurs utilisés (indiqués dans le sous-titre de chaque
graphique et tableau), en particulier :

Makowski D, Ben-Shachar MS, Patil I, Lüdecke D (2020). Methods and
algorithms for correlation analysis in R. *Journal of Open Source
Software* 5(51):2306.
[doi:10.21105/joss.02306](https://doi.org/10.21105/joss.02306)

## Licence

MIT © Clément Poupelin ([ORCID
0009-0005-7308-0368](https://orcid.org/0009-0005-7308-0368))
