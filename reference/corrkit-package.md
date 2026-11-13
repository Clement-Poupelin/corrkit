# corrkit: correlation analyses with annotated heatmaps and tables

corrkit does not compute correlations itself: it calls validated engines
and adds a way of working on top of them.

## Details

**Engines** (an existing function for every computation, with a base-R
fallback when a suggested package is missing)

- one pair, unadjusted:
  [`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html)

- one pair, partial:
  [`ppcor::pcor.test()`](https://rdrr.io/pkg/ppcor/man/pcor.test.html)
  (Kim 2015, doi:10.5351/CSAM.2015.22.6.665)

- matrix of all pairs:
  [`psych::corr.test()`](https://rdrr.io/pkg/psych/man/corr.test.html)
  or
  [`correlation::correlation()`](https://easystats.github.io/correlation/reference/correlation.html)
  (Makowski et al. 2020, doi:10.21105/joss.02306), always with raw
  p-values; the adjustment is done by corrkit with
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) over an
  explicit family

- binary effect sizes:
  [`stats::t.test()`](https://rdrr.io/r/stats/t.test.html),
  [`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html),
  [`effectsize::rank_biserial()`](https://easystats.github.io/effectsize/reference/rank_biserial.html)
  for the confidence interval

**What corrkit adds**

- an explicit multiple-testing family (`adjust_by`, BH by default);

- a dedicated treatment of binary (0/1) variables: separate blocks in
  the plots, a `[0/1]` marker, a reading note
  ([`cor_binary_note()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_note.md))
  and the effect sizes to report instead of the correlation
  ([`cor_binary_effects()`](https://clement-poupelin.github.io/corrkit/reference/cor_binary_effects.md));

- outputs whose subtitle, caption or source notes state the method, the
  engine, the formulas, the adjustment and how to read them:
  [`plot_cor_heatmap()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_heatmap.md),
  [`plot_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_cross.md),
  [`plot_cor_pair()`](https://clement-poupelin.github.io/corrkit/reference/plot_cor_pair.md),
  [`tbl_cor_matrix()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_matrix.md),
  [`tbl_cor_cross()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_cross.md),
  [`tbl_cor_long()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_long.md),
  [`tbl_cor_binary()`](https://clement-poupelin.github.io/corrkit/reference/tbl_cor_binary.md).

Start with
[`vignette("corrkit")`](https://clement-poupelin.github.io/corrkit/articles/corrkit.md).

## See also

Useful links:

- <https://clement-poupelin.github.io/corrkit/>

- <https://github.com/clement-poupelin/corrkit>

- Report bugs at <https://github.com/clement-poupelin/corrkit/issues>

## Author

**Maintainer**: Clément Poupelin <clementjc.poupelin@gmail.com>
([ORCID](https://orcid.org/0009-0005-7308-0368))

Authors:

- Clément Poupelin <clementjc.poupelin@gmail.com>
  ([ORCID](https://orcid.org/0009-0005-7308-0368))
