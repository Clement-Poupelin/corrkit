# ====================================================================
# notes.R -- reading notes, captions and source notes
# ====================================================================
# Every figure and every table states, in this order:
#   (1) the p-values: raw or adjusted, method, family of the adjustment;
#   (2) the coefficient: formula (notation of Makowski et al. 2020), test, df;
#   (3) the confidence interval (tables);
#   (4) binary variables: coefficient actually computed, how to read it;
#   (5) the engine(s): function() {package}.
# The same helpers feed the ggplot2 captions and the gt source notes, so a
# figure and its table never disagree.
# ====================================================================

#' Palette of the correlation heatmaps and tables (named low / high, white at 0)
#' @export
COR_PAL <- c(low = "#2166AC", high = "#B2182B")

#' The three methods with a full treatment (partial correlation, notes, formulas)
#' @noRd
.core_methods <- c("spearman", "pearson", "kendall")

#' Methods offered by correlation::correlation() that corrkit passes through (matrix only)
#' @noRd
.extra_methods <- c("biweight", "percentage", "blomqvist", "hoeffding", "gaussian", "shepherd", "distance")


# --------------------------------------------------------------------
# Reading note for binary variables
# --------------------------------------------------------------------

#' Reading note for a correlation involving a binary (0/1) variable
#'
#' * **Pearson** + binary = point-biserial r_pb. Its test is Student's t-test
#'   (equal variances); partial (adjusted) = slope test of
#'   `lm(Y ~ X + covariates)` (Frisch-Waugh-Lovell: same t, same df). |r_pb|
#'   depends on the split p/q: for a biomarker categorised from a continuous
#'   measurement, r_pb = r x phi(c) / sqrt(p q) (c = cut-off in SD units), i.e.
#'   at most 0.80 x r at the median, 0.70 x r for a 20/80 split.
#' * **Spearman** + binary = Pearson on ranks; a linear function of the
#'   Mann-Whitney U (same test, approximately). r_s = r_rb x sqrt(3 p q) (no
#'   ties in Y), with r_rb = 2 AUC - 1 the rank-biserial correlation:
#'   |r_s| <= 0.87 (50/50), 0.69 (20/80).
#' * **Kendall** + binary: the numerator of tau (concordant - discordant pairs)
#'   is n_0 n_1 r_rb, i.e. the Mann-Whitney comparison; its maximum depends on
#'   the split and on the ties.
#'
#' In every case the effect size to report is r_rb or the difference between
#' groups ([cor_binary_effects()]), not the correlation.
#'
#' @param meth `"spearman"`, `"pearson"` or `"kendall"`
#' @param lang `"en"` (figure captions, table notes) or `"fr"` (callouts of a French report)
#' @param short `TRUE` = one line for a caption; `FALSE` = full note with formulas and references
#' @param mark marker of binary variables used in the text (same as `bin_mark` of the plots)
#' @return a character string
#' @references
#' Kerby DS (2014). The simple difference formula: an approach to teaching
#' nonparametric correlation. *Comprehensive Psychology* 3:11.IT.3.1. doi:10.2466/11.IT.3.1
#'
#' Babchishin KM, Helmus LM (2016). The influence of base rates on
#' correlations. *Behavior Research Methods* 48:1021-1031. doi:10.3758/s13428-015-0627-7
#' @export
#' @examples
#' cor_binary_note("spearman", short = TRUE)
#' cat(cor_binary_note("pearson", lang = "fr"))
cor_binary_note <- function(meth = c("spearman", "pearson", "kendall"), lang = c("en", "fr"), short = FALSE,
                            mark = " [0/1]") {
  meth <- match.arg(meth); lang <- match.arg(lang)
  gsub("{MARK}", trimws(mark), .cor_binary_text(meth, lang, short), fixed = TRUE)
}

# Internal: texts of cor_binary_note(), with the placeholder {MARK} for the marker
.cor_binary_text <- function(meth, lang, short) {
  refs <- paste("Kerby 2014, Comprehensive Psychology, doi:10.2466/11.IT.3.1;",
                "Babchishin & Helmus 2016, Behav Res Methods, doi:10.3758/s13428-015-0627-7")
  en <- lang == "en"
  if (short) {
    return(switch(meth,
      pearson = if (en)
        "{MARK} binary: point-biserial r (= t-test / lm slope); |r| shrinks with an unbalanced split, not comparable with continuous variables"
        else "{MARK} binaire : r point-bis\u00e9rial (= test t / pente du lm) ; |r| diminue si les groupes sont d\u00e9s\u00e9quilibr\u00e9s, non comparable aux variables continues",
      spearman = if (en)
        "{MARK} binary: r_s \u2248 Mann-Whitney (rank difference between groups); |r_s| \u2264 \u221a(3pq) \u2264 0.87, not comparable with continuous variables"
        else "{MARK} binaire : r_s \u2248 Mann-Whitney (\u00e9cart de rangs entre groupes) ; |r_s| \u2264 \u221a(3pq) \u2264 0,87, non comparable aux variables continues",
      kendall = if (en)
        "{MARK} binary: \u03c4 \u2248 Mann-Whitney (concordant \u2212 discordant pairs between groups); bounded by the split, not comparable with continuous variables"
        else "{MARK} binaire : \u03c4 \u2248 Mann-Whitney (paires concordantes \u2212 discordantes entre groupes) ; born\u00e9 par la r\u00e9partition, non comparable aux variables continues"))
  }
  core <- switch(meth,
    pearson = if (en) paste(
      "Binary variables (marked {MARK}, upper class = 1): Pearson r is the point-biserial correlation.",
      "Its test is Student's t-test between the two groups (equal variances) and, adjusted for covariates,",
      "the test of the slope of lm(score ~ group + covariates): same t, same degrees of freedom.",
      "Its size depends on the split p/q between the groups: for a biomarker categorised from a continuous",
      "measurement, r_pb = r \u00d7 \u03c6(c) / \u221a(pq) (c = cut-off in SD units), i.e. at most 0.80 \u00d7 r at the median",
      "and 0.70 \u00d7 r for a 20/80 split. It is therefore not comparable with the r of a continuous biomarker:",
      "read the size of the effect as the difference in mean score between groups.")
    else paste(
      "Variables binaires (not\u00e9es {MARK}, classe haute = 1) : le r de Pearson est la corr\u00e9lation point-bis\u00e9riale.",
      "Son test est le test t de Student entre les deux groupes (variances \u00e9gales) et, ajust\u00e9 sur des covariables,",
      "le test de la pente de lm(score ~ groupe + covariables) : m\u00eame t, m\u00eames degr\u00e9s de libert\u00e9.",
      "Sa valeur d\u00e9pend de la r\u00e9partition p/q entre groupes : pour un biomarqueur cat\u00e9goris\u00e9 \u00e0 partir d'une mesure continue,",
      "r_pb = r \u00d7 \u03c6(c) / \u221a(pq) (c = seuil en unit\u00e9s d'\u00e9cart-type), soit au plus 0,80 \u00d7 r \u00e0 la m\u00e9diane et",
      "0,70 \u00d7 r pour une r\u00e9partition 20/80. Il n'est donc pas comparable au r d'un biomarqueur continu :",
      "la taille de l'effet se lit sur la diff\u00e9rence de score moyen entre groupes."),
    spearman = if (en) paste(
      "Binary variables (marked {MARK}, upper class = 1): Spearman r_s is a Pearson correlation on ranks;",
      "with two groups it is a linear function of the Mann-Whitney U (same test, approximately), i.e. it measures",
      "whether scores tend to rank higher in the upper group. r_s = r_rb \u00d7 \u221a(3pq), where r_rb = 2 AUC \u2212 1 is the",
      "rank-biserial correlation (p/q = split between groups, no ties): |r_s| \u2264 0.87 at 50/50 and 0.69 at 20/80.",
      "It is therefore not comparable with the r_s of a continuous biomarker: report r_rb or the difference",
      "between groups as the effect size. The partial (adjusted) version is read the same way.")
    else paste(
      "Variables binaires (not\u00e9es {MARK}, classe haute = 1) : le r_s de Spearman est une corr\u00e9lation de Pearson sur les rangs ;",
      "avec deux groupes, c'est une fonction lin\u00e9aire du U de Mann-Whitney (m\u00eame test, \u00e0 l'approximation pr\u00e8s) : il mesure",
      "si les scores tendent \u00e0 \u00eatre mieux class\u00e9s dans le groupe haut. r_s = r_rb \u00d7 \u221a(3pq), o\u00f9 r_rb = 2 AUC \u2212 1 est la",
      "corr\u00e9lation rang-bis\u00e9riale (p/q = r\u00e9partition entre groupes, sans ex aequo) : |r_s| \u2264 0,87 \u00e0 50/50 et 0,69 \u00e0 20/80.",
      "Il n'est donc pas comparable au r_s d'un biomarqueur continu : la taille de l'effet est r_rb ou la diff\u00e9rence",
      "entre groupes. La version partielle (ajust\u00e9e) se lit de la m\u00eame fa\u00e7on."),
    kendall = if (en) paste(
      "Binary variables (marked {MARK}, upper class = 1): the numerator of Kendall \u03c4 (concordant \u2212 discordant pairs)",
      "only counts pairs taken in different groups and equals n_0 \u00d7 n_1 \u00d7 r_rb, where r_rb = 2 AUC \u2212 1 is the",
      "rank-biserial correlation: it is the Mann-Whitney comparison. Its maximum depends on the split p/q and on the",
      "ties, so it is not comparable with the \u03c4 of a continuous biomarker: report r_rb or the difference between",
      "groups as the effect size.")
    else paste(
      "Variables binaires (not\u00e9es {MARK}, classe haute = 1) : le num\u00e9rateur du \u03c4 de Kendall (paires concordantes \u2212 discordantes)",
      "ne compte que les paires prises dans des groupes diff\u00e9rents et vaut n_0 \u00d7 n_1 \u00d7 r_rb, o\u00f9 r_rb = 2 AUC \u2212 1 est la",
      "corr\u00e9lation rang-bis\u00e9riale : c'est la comparaison de Mann-Whitney. Son maximum d\u00e9pend de la r\u00e9partition p/q et des",
      "ex aequo : il n'est donc pas comparable au \u03c4 d'un biomarqueur continu ; la taille de l'effet est r_rb ou la",
      "diff\u00e9rence entre groupes."))
  paste0(core, " (", refs, ").")
}


# --------------------------------------------------------------------
# Building blocks of captions and source notes
# --------------------------------------------------------------------

# Internal: caption made of blocks (list of character vectors), one blank line between blocks;
# each line wrapped at `width` characters (NULL = no wrapping)
.wrap_caption <- function(blocks, width = NULL) {
  if (!is.list(blocks)) blocks <- list(blocks)
  wrap1 <- function(l) if (is.null(width)) l else paste(strwrap(l, width), collapse = "\n")
  txt <- vapply(blocks, function(b) {
    b <- b[!is.na(b) & nzchar(b)]
    if (length(b)) paste(vapply(b, wrap1, character(1)), collapse = "\n") else ""
  }, character(1))
  paste(txt[nzchar(txt)], collapse = "\n\n")
}

# Internal: full name of a p.adjust() method
.p_adjust_name <- function(m) switch(m, BH = , fdr = "Benjamini-Hochberg (FDR)",
                                     BY = "Benjamini-Yekutieli (FDR)", holm = "Holm (FWER)",
                                     hochberg = "Hochberg (FWER)", hommel = "Hommel (FWER)",
                                     bonferroni = "Bonferroni (FWER)", none = "none", m)

# Internal: p-value line
.cor_p_line <- function(p_adjust, family_txt) {
  if (identical(p_adjust, "none") || is.null(family_txt))
    return("Raw p (no multiple-testing adjustment): * < 0.05, ** < 0.01, *** < 0.001; no star = p \u2265 0.05.")
  sprintf("Adjusted p, %s over %s: * < 0.05, ** < 0.01, *** < 0.001; no star = adjusted p \u2265 0.05.",
          .p_adjust_name(p_adjust), family_txt)
}

# Internal: symbol of the coefficient, plotmath (legend) and text
.cor_sym <- function(meth) switch(meth, spearman = quote(italic(r)[s]), kendall = quote(tau), quote(italic(r)))
.cor_sym_txt <- function(meth) switch(meth, spearman = "r_s", kendall = "\u03c4", "r")
# Internal: symbol for gt column labels (Unicode subscript s)
.cor_sym_tbl <- function(meth) switch(meth, spearman = "r\u209b", kendall = "\u03c4", "r")

# Internal: title-case name of the method
.meth_name <- function(meth) switch(meth, spearman = "Spearman", pearson = "Pearson", kendall = "Kendall",
                                    percentage = "Percentage bend", biweight = "Biweight midcorrelation",
                                    blomqvist = "Blomqvist", hoeffding = "Hoeffding", gaussian = "Gaussian rank",
                                    shepherd = "Shepherd's pi", distance = "Distance", tools::toTitleCase(meth))

# Internal: formula line(s), notation of Makowski et al. 2020
# (r_xy Pearson, r_s Spearman, tau Kendall, r_xy.z partial with e_x.z = residuals of x on z)
.cor_formula_lines <- function(meth, partial = FALSE, k = 0) {
  df <- if (partial) "n \u2212 2 \u2212 k" else "n \u2212 2"
  if (meth == "spearman") {
    coef <- "Continuous: r_s = cov(rank_x, rank_y) / (SD(rank_x) \u00d7 SD(rank_y))"
    if (partial) coef <- paste0(coef, "; partial: r_s,xy.z = r(e_x.z, e_y.z), residuals of the ranks on the k = ", k, " covariate terms z")
    return(paste0(coef, "; t-test, ", df, " df."))
  }
  if (meth == "pearson") {
    coef <- "Continuous: r_xy = cov(x, y) / (SD_x \u00d7 SD_y)"
    if (partial) coef <- paste0(coef, "; partial: r_xy.z = r(e_x.z, e_y.z), residuals on the k = ", k, " covariate terms z")
    return(paste0(coef, "; t-test, ", df, " df."))
  }
  if (meth == "kendall") {
    coef <- "Continuous: \u03c4_b = (concordant \u2212 discordant pairs) / \u221a(pairs untied on x \u00d7 pairs untied on y)"
    if (partial) coef <- paste0(coef, "; partial: \u03c4_xy.z from the inverse of the \u03c4 matrix of (x, y, z), k = ", k, " covariate terms z")
    return(paste0(coef, "; z-test (normal approximation)."))
  }
  sprintf("Coefficient: %s correlation as defined in correlation::correlation() (Makowski et al. 2020); test as returned by that function.",
          .meth_name(meth))
}

# Internal: binary-variable line (formula of the coefficient actually computed, and reading)
# Pearson : point-biserial r_pb = (M_1 - M_0) / SD_y x sqrt(pq) (SD with denominator n)
# Spearman: r_s = r_rb x sqrt(3pq), r_rb = 2 AUC - 1 rank-biserial correlation (Mann-Whitney), no ties
.cor_binary_lines <- function(meth, bin_mark = " [0/1]") {
  m <- trimws(bin_mark)
  switch(meth,
    pearson = sprintf("Binary %s: point-biserial r_pb = (M_1 \u2212 M_0) / SD_y \u00d7 \u221a(pq), M = group means, p, q = group proportions (\u2261 Student's t-test); depends on p, q, not comparable with continuous pairs.", m),
    spearman = sprintf("Binary %s: r_s = r_rb \u00d7 \u221a(3pq), r_rb = 2 AUC \u2212 1 (rank-biserial, \u2261 Mann-Whitney test), p, q = group proportions; |r_s| \u2264 0.87, not comparable with continuous pairs.", m),
    kendall = sprintf("Binary %s: numerator of \u03c4 = n_0 n_1 \u00d7 r_rb, r_rb = 2 AUC \u2212 1 (rank-biserial, \u2261 Mann-Whitney test); bounded by the split p, q, not comparable with continuous pairs.", m),
    sprintf("Binary %s: coefficient computed on a 0/1 variable; depends on the split between groups, not comparable with continuous pairs.", m))
}

# Internal: confidence-interval line (tables)
.cor_ci_line <- function(ci, meth, ci_engine = NULL, partial = FALSE) {
  if (is.null(ci) || is.na(ci)) return(NULL)
  lvl <- sprintf("%g%% CI", 100 * ci)
  if (!meth %in% .core_methods) return(paste0(lvl, ": as returned by correlation::correlation() {correlation}."))
  nk <- if (partial) " \u2212 k" else ""
  se <- switch(meth, spearman = sprintf("SE = \u221a(1.06 / (n \u2212 3%s)) (Fieller et al. 1957)", nk),
               kendall = sprintf("SE = \u221a(0.437 / (n \u2212 4%s)) (Fieller et al. 1957)", nk),
               sprintf("SE = 1 / \u221a(n \u2212 3%s)", nk))
  eng <- if (identical(ci_engine, "correlation::cor_to_ci")) "correlation::cor_to_ci() {correlation}" else "base R"
  sprintf("%s: Fisher z transformation, %s; computed with %s.", lvl, se, eng)
}

# Internal: "function() {package}" text of the engines used
.cor_engine_txt <- function(e) {
  e <- unique(e[!is.na(e) & e != "none"])
  if (!length(e)) return("no engine (no usable pair)")
  txt <- vapply(e, function(x) switch(x,
    "stats::cor.test"         = "stats::cor.test() {stats}",
    "ppcor::pcor.test"        = "ppcor::pcor.test() {ppcor}",
    "psych::corr.test"        = "psych::corr.test() {psych}",
    "correlation::correlation" = "correlation::correlation() {correlation}",
    "residuals (base R)"      = "residuals of stats::lm.fit() + stats::cor() {stats}", x), character(1))
  paste(txt, collapse = " + ")
}

# Internal: notation / reference line of the tables
.cor_ref_line <- function()
  "Notation: Makowski D, Ben-Shachar MS, Patil I, L\u00fcdecke D (2020), J Open Source Softw 5(51):2306, doi:10.21105/joss.02306."

# Internal: more room around the subtitle and the caption; with ggtext (see
# .use_markdown()), subtitle, caption and axis titles are markdown elements so
# that the notation prepared by .labs_safe() is typeset (subscripts, italics)
.cor_text_theme <- function() {
  if (.use_markdown())
    return(ggplot2::theme(
      plot.subtitle = ggtext::element_markdown(hjust = 0, halign = 0, lineheight = 1.15),
      plot.caption  = ggtext::element_markdown(hjust = 0, halign = 0, lineheight = 1.2,
                                               margin = ggplot2::margin(t = 14)),
      axis.title.x  = ggtext::element_markdown(),
      axis.title.y  = ggtext::element_markdown()))
  ggplot2::theme(plot.subtitle = ggplot2::element_text(lineheight = 1.15),
                 plot.caption = ggplot2::element_text(hjust = 0, lineheight = 1.2,
                                                      margin = ggplot2::margin(t = 14)))
}
