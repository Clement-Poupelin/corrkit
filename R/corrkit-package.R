#' corrkit: correlation analyses with annotated heatmaps and tables
#'
#' corrkit does not compute correlations itself: it calls validated engines and
#' adds a way of working on top of them.
#'
#' **Engines** (an existing function for every computation, with a base-R
#' fallback when a suggested package is missing)
#' * one pair, unadjusted: [stats::cor.test()]
#' * one pair, partial: `ppcor::pcor.test()` (Kim 2015, doi:10.5351/CSAM.2015.22.6.665)
#' * matrix of all pairs: `psych::corr.test()` or `correlation::correlation()`
#'   (Makowski et al. 2020, doi:10.21105/joss.02306), always with raw p-values;
#'   an optional adjustment is done by corrkit with [stats::p.adjust()] over an explicit family
#' * binary effect sizes: [stats::t.test()], [stats::wilcox.test()],
#'   `effectsize::rank_biserial()` for the confidence interval
#'
#' **What corrkit adds**
#' * raw p-values by default and, on request, an adjustment over an explicit
#'   family of tests (`p_adjust`, `adjust_by`);
#' * a dedicated treatment of binary (0/1) variables: separate blocks in the
#'   plots, a `[0/1]` marker, a reading note ([cor_binary_note()]) and the
#'   effect sizes to report instead of the correlation ([cor_binary_effects()]);
#' * outputs whose subtitle, caption or source notes state the method, the
#'   engine, the formulas, the adjustment and how to read them:
#'   [plot_cor_heatmap()], [plot_cor_cross()], [plot_cor_pair()],
#'   [tbl_cor_matrix()], [tbl_cor_cross()], [tbl_cor_long()], [tbl_cor_binary()].
#'
#' Start with `vignette("corrkit")`.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom dplyr %>%
#' @importFrom utils globalVariables
NULL

# Column names used inside dplyr / ggplot2 calls (avoid "no visible binding" notes)
utils::globalVariables(c(
  "Var1", "Var2", "X", "Y", "res", "x_type", "y_type", "rho", "p", "p_adj",
  "x", "y", "lab", "Xl", "Yl", "i", "j", "x_block", "y_block", "block",
  "ex", "ey", "g", "v", "ci_low", "ci_high", "n"
))
