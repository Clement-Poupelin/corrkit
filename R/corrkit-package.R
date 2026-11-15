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
#'   [plot_cor_heatmap()], [plot_cor_forest()], [plot_cor_pair()] (or `plot()`),
#'   [tbl_cor()], [tbl_cor_binary()].
#'
#' **One computing function**: [cor_compute()] (matrix, cross or single pair,
#' simple or partial), whose result has `print()`, `summary()`, `plot()`,
#' `autoplot()`, `as.matrix()` and `as.data.frame()` methods.
#'
#' **Typesetting of the plot texts** (options)
#' * `options(corrkit.markdown = TRUE)` (default): when the ggtext package is
#'   installed, subtitles, captions and axis titles are rendered as markdown,
#'   so that the notation is typeset as in LaTeX (r_s becomes an italic r with
#'   a subscript s, r_xy.z, e_x.z, r_rb, M_1, SD_y likewise). `FALSE` gives plain text.
#' * `options(corrkit.unicode = FALSE)` (default): the square root, the minus
#'   sign, tau, <= ... are written `sqrt`, `-`, `tau`, `<=`, because the pdf()
#'   device (used by R CMD check and by `ggsave("figure.pdf")`) cannot draw
#'   them. Set it to `TRUE` if you save your figures as PNG, or as PDF with
#'   `ggsave(..., device = cairo_pdf)`. The legend of the colour scale is
#'   always typeset (plotmath), on every device.
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
  "ex", "ey", "g", "v", "ci_low", "ci_high", "n", "r", "key", "lo", "hi", "sig", "dir", "panel"
))
