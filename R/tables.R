# ====================================================================
# tables.R -- gt tables with informative subtitles and source notes
# ====================================================================
# Every table carries:
#   title / subtitle : method, adjustment (covariates), N, engine(s)
#   cells            : coefficient + stars (colour = same scale as the heatmaps)
#   footnotes        : on the binary [0/1] variables
#   source notes     : p-values, formula, CI, binary reading, notation reference
# ====================================================================

# Internal: add the non-empty lines as source notes
.gt_notes <- function(g, notes) {
  notes <- unlist(notes)
  notes <- notes[!is.na(notes) & nzchar(notes)]
  for (nt in notes) g <- gt::tab_source_note(g, source_note = nt)
  g
}

# Internal: colour the text column `target` from the numeric column `source` (scale -1..1)
.gt_colour <- function(g, source, target, palette = COR_PAL) {
  for (i in seq_along(source))
    g <- gt::data_color(g, columns = source[i], target_columns = target[i], method = "numeric",
                        palette = c(palette[["low"]], "#FFFFFF", palette[["high"]]), domain = c(-1, 1),
                        na_color = "#FFFFFF")
  g
}

# Internal: common style
.gt_style <- function(g)
  gt::tab_options(g, table.font.size = gt::px(13), heading.align = "left",
                  source_notes.font.size = gt::px(11), footnotes.font.size = gt::px(11),
                  heading.subtitle.font.size = gt::px(12))


#' Table of cross correlations X x Y (wide: one row per predictor, one column per score)
#'
#' Same content and same colour scale as [plot_cor_cross()]. Binary
#' predictors are grouped under "Binary [0/1]" and carry a footnote.
#'
#' @param x output of [cor_cross()] or [plot_cor_cross()]
#' @param labels named vector variable -> label (default: those given to [cor_cross()])
#' @param digits decimals of the coefficients
#' @param palette colours (named low / high); `NULL` = no colour
#' @param show_note add the reading lines of binary variables
#' @param title,subtitle override the automatic texts
#' @param stub_label header of the row labels
#' @return a `gt_tbl`
#' @export
#' @examples
#' d  <- simulate_cor_data(150)
#' cc <- cor_cross(d, c("CRP", "IL6", "BIN_CRP"), c("depression", "fatigue"),
#'                 covariates = "age", labels = attr(d, "labels"))
#' tbl_cor_cross(cc)
tbl_cor_cross <- function(x, labels = NULL, digits = 2, palette = COR_PAL, show_note = TRUE,
                          title = NULL, subtitle = NULL, stub_label = "Predictor") {
  cc <- .as_cor_cross(x); tab <- cc$table; m <- cc$meta
  labels <- labels %||% m$labels
  bin_mark <- m$bin_mark %||% " [0/1]"
  x_vars <- m$x_vars; y_vars <- m$y_vars
  x_lab <- mark_binary(get_label(x_vars, labels), m$x_types, bin_mark)
  y_lab <- mark_binary(get_label(y_vars, labels), m$y_types, bin_mark)
  fmt <- paste0("%.", digits, "f")
  tab$lab <- ifelse(is.na(tab$rho), "", paste0(sprintf(fmt, tab$rho), p_to_stars(tab$p_adj, ns = "")))

  lab_cols <- paste0("lab__", seq_along(y_vars)); r_cols <- paste0("r__", seq_along(y_vars))
  wide <- data.frame(row_lab = x_lab, stringsAsFactors = FALSE)
  for (j in seq_along(y_vars)) {
    sub <- tab[tab$Y == y_vars[j], ]
    wide[[lab_cols[j]]] <- sub$lab[match(x_vars, sub$X)]
    wide[[r_cols[j]]]   <- sub$rho[match(x_vars, sub$X)]
  }
  is_bin <- m$x_types == "binary"
  bin_block <- paste0("Binary", bin_mark)
  wide$group <- ifelse(is_bin, bin_block, "Continuous")
  split_rows <- any(is_bin) && any(!is_bin)
  ord <- order(is_bin)                      # continuous first, original order kept within blocks
  wide <- wide[ord, , drop = FALSE]; is_bin <- is_bin[ord]

  g <- if (split_rows) gt::gt(wide, rowname_col = "row_lab", groupname_col = "group")
       else gt::gt(wide, rowname_col = "row_lab")
  if (!split_rows) g <- gt::cols_hide(g, "group")
  if (!is.null(palette)) g <- .gt_colour(g, r_cols, lab_cols, palette)
  g <- g %>%
    gt::cols_hide(r_cols) %>%
    gt::cols_label(.list = stats::setNames(as.list(y_lab), lab_cols)) %>%
    gt::cols_align("center", columns = lab_cols) %>%
    gt::tab_stubhead(stub_label) %>%
    gt::tab_header(title = title %||% sprintf("%s correlations", .meth_name(m$meth)),
                   subtitle = subtitle %||% .cor_cross_subtitle(m, tab, labels, newline = " | "))
  y_bin <- m$y_types == "binary"
  if (any(y_bin) && any(!y_bin))
    g <- g %>% gt::tab_spanner("Continuous", columns = lab_cols[!y_bin], id = "sp_cont") %>%
      gt::tab_spanner(bin_block, columns = lab_cols[y_bin], id = "sp_bin")
  if (any(is_bin) && show_note)
    g <- gt::tab_footnote(g, footnote = cor_binary_note(m$meth, "en", short = TRUE, mark = bin_mark),
                          locations = gt::cells_stub(rows = which(is_bin)))
  g <- .gt_notes(g, list(
    paste0("Cells: ", .cor_sym_txt(m$meth), " and stars; colour = same scale as the heatmap (blue < 0 < red). ",
           "Confidence intervals: tbl_cor_long()."),
    .cor_p_line(m$p_adjust, m$family),
    .cor_formula_lines(m$meth, m$partial, m$k),
    if (show_note && any(c(m$x_types, m$y_types) == "binary")) .cor_binary_lines(m$meth, bin_mark),
    .cor_ref_line()))
  .gt_style(g)
}


#' Long table of correlations: one row per pair, with CI, raw and adjusted p
#'
#' Works for [cor_cross()] / [plot_cor_cross()] outputs (rows grouped by the
#' family of the adjustment) and for [cor_test_matrix()] outputs (all pairs).
#' Rows with a p-value (adjusted, if an adjustment was requested) below
#' `alpha` are in bold; the adjusted column is hidden when p-values are raw.
#'
#' @param x output of [cor_cross()], [plot_cor_cross()] or [cor_test_matrix()]
#' @param labels named vector variable -> label
#' @param digits decimals of the coefficients
#' @param alpha threshold for bold rows
#' @param palette colours (named low / high) of the coefficient column; `NULL` = no colour
#' @param title,subtitle override the automatic texts
#' @param sort_by `"none"` (order of the variables) or `"abs"` (decreasing |coefficient| within groups)
#' @return a `gt_tbl`
#' @export
#' @examples
#' d  <- simulate_cor_data(150)
#' cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "BIN_CRP"))
#' tbl_cor_long(cm, labels = attr(d, "labels"))
tbl_cor_long <- function(x, labels = NULL, digits = 2, alpha = 0.05, palette = COR_PAL,
                         title = NULL, subtitle = NULL, sort_by = c("none", "abs")) {
  sort_by <- match.arg(sort_by)
  bin_mark <- " [0/1]"
  if (inherits(x, "cor_matrix")) {
    meth <- x$method; L <- x$long
    d <- data.frame(var1 = L$Var1, var2 = L$Var2, type1 = L$type1, type2 = L$type2, n = L$n,
                    r = L$r, ci_low = L$ci_low, ci_high = L$ci_high, p = L$p, p_adj = L$p_adj,
                    group = "", stringsAsFactors = FALSE)
    head1 <- "Variable 1"; head2 <- "Variable 2"
    p_line <- .cor_p_line(x$p_adjust, sprintf("the %d distinct pairs", nrow(L)))
    formula <- .cor_formula_lines(meth)
    ci_line <- if (identical(x$ci_engine, "correlation::correlation"))
      sprintf("%g%% CI: as returned by correlation::correlation() {correlation}.", 100 * (x$ci %||% 0.95))
      else .cor_ci_line(x$ci, meth, x$ci_engine)
    sub_default <- sprintf("%s correlation, all pairs | pairwise complete N = %s | Computed with %s",
                           .meth_name(meth), .n_range_txt(L$n), .cor_engine_txt(x$engine))
    title_default <- sprintf("%s correlations", .meth_name(meth))
    has_bin <- any(x$types == "binary")
  } else {
    cc <- .as_cor_cross(x); tab <- cc$table; m <- cc$meta
    labels <- labels %||% m$labels; bin_mark <- m$bin_mark %||% bin_mark
    meth <- m$meth
    d <- data.frame(var1 = tab$X, var2 = tab$Y, type1 = tab$x_type, type2 = tab$y_type, n = tab$n,
                    r = tab$rho, ci_low = tab$ci_low, ci_high = tab$ci_high, p = tab$p, p_adj = tab$p_adj,
                    stringsAsFactors = FALSE)
    d$group <- switch(m$adjust_by, y = get_label(tab$Y, labels), x = get_label(tab$X, labels), "")
    head1 <- "Predictor"; head2 <- "Score"
    p_line <- .cor_p_line(m$p_adjust, m$family)
    formula <- .cor_formula_lines(meth, m$partial, m$k)
    ci_line <- .cor_ci_line(m$ci, meth, m$ci_engine, m$partial)
    sub_default <- .cor_cross_subtitle(m, tab, labels, newline = " | ")
    title_default <- sprintf("%s correlations", .meth_name(meth))
    has_bin <- any(c(m$x_types, m$y_types) == "binary")
  }
  if (sort_by == "abs") d <- d[order(match(d$group, unique(d$group)), -abs(d$r)), , drop = FALSE]
  fmt <- paste0("%.", digits, "f")
  short_type <- function(t) ifelse(t == "binary", "bin", "cont")
  d$lab1 <- mark_binary(get_label(d$var1, labels), d$type1, bin_mark)
  d$lab2 <- mark_binary(get_label(d$var2, labels), d$type2, bin_mark)
  d$pair <- paste(short_type(d$type1), "\u00d7", short_type(d$type2))
  d$coef <- ifelse(is.na(d$r), "", sprintf(fmt, d$r))
  d$ci   <- ifelse(is.na(d$ci_low), "", sprintf(paste0("[", fmt, ", ", fmt, "]"), d$ci_low, d$ci_high))
  d$p_txt    <- .fmt_p(d$p)
  raw_only <- identical(if (inherits(x, "cor_matrix")) x$p_adjust else m$p_adjust, "none")
  if (raw_only) d$p_txt <- paste0(d$p_txt, ifelse(is.na(d$p), "", paste0(" ", p_to_stars(d$p, ns = ""))))
  d$padj_txt <- paste0(.fmt_p(d$p_adj), ifelse(is.na(d$p_adj), "", paste0(" ", p_to_stars(d$p_adj, ns = ""))))
  grouped <- any(nzchar(d$group))
  # the variable that defines the row groups is not repeated in a column
  drop_col <- if (!grouped) NULL else if (identical(m$adjust_by, "x")) "lab1" else "lab2"
  show <- data.frame(group = d$group, lab1 = d$lab1, lab2 = d$lab2, pair = d$pair, n = d$n, coef = d$coef,
                     r_num = d$r, ci = d$ci, p_txt = d$p_txt, padj_txt = d$padj_txt, stringsAsFactors = FALSE)
  if (!is.null(drop_col)) show[[drop_col]] <- NULL
  g <- if (grouped) gt::gt(show, groupname_col = "group") else gt::gt(show)
  if (!grouped) g <- gt::cols_hide(g, "group")
  if (!is.null(palette)) g <- .gt_colour(g, "r_num", "coef", palette)
  ci_lab <- if (any(nzchar(d$ci))) sprintf("%g%% CI", 100 * ((if (inherits(x, "cor_matrix")) x$ci else m$ci) %||% 0.95)) else "CI"
  labs <- list(pair = "Pair", n = "n", coef = .cor_sym_tbl(meth), ci = ci_lab,
               p_txt = "p (raw)", padj_txt = "p (adjusted)")
  if (!identical(drop_col, "lab1")) labs$lab1 <- head1
  if (!identical(drop_col, "lab2")) labs$lab2 <- head2
  g <- g %>%
    gt::cols_hide("r_num") %>%
    gt::cols_label(.list = labs) %>%
    gt::cols_align("center", columns = c("pair", "n", "coef", "ci", "p_txt", "padj_txt")) %>%
    gt::tab_header(title = title %||% title_default, subtitle = subtitle %||% sub_default)
  sig <- which(!is.na(d$p_adj) & d$p_adj < alpha)
  if (length(sig)) g <- gt::tab_style(g, style = gt::cell_text(weight = "bold"), locations = gt::cells_body(rows = sig))
  if (!any(nzchar(d$ci))) g <- gt::cols_hide(g, "ci")
  if (raw_only) g <- gt::cols_hide(g, "padj_txt")
  g <- .gt_notes(g, list(
    if (grouped) sprintf("Row groups = family of the p adjustment (%s).", if (!inherits(x, "cor_matrix")) m$family else ""),
    sprintf("Pair: cont = continuous, bin = binary%s. Bold = %s < %g.", trimws(bin_mark),
            if (raw_only) "p" else "adjusted p", alpha),
    p_line, formula, ci_line,
    if (has_bin) .cor_binary_lines(meth, bin_mark),
    .cor_ref_line()))
  .gt_style(g)
}


#' Table of a correlation matrix (lower triangle), same layout as the heatmap
#'
#' @param cm output of [cor_test_matrix()]
#' @param labels named vector variable -> label
#' @param use_adj stars on adjusted (`TRUE`) or raw p-values
#' @param order `"original"` or `"hclust"` (same ordering as [plot_cor_heatmap()])
#' @param digits decimals
#' @param palette colours (named low / high); `NULL` = no colour
#' @param bin_mark marker of binary variables
#' @param title,subtitle override the automatic texts
#' @return a `gt_tbl`
#' @export
#' @examples
#' d  <- simulate_cor_data(150)
#' cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"))
#' tbl_cor_matrix(cm, labels = attr(d, "labels"))
tbl_cor_matrix <- function(cm, labels = NULL, use_adj = TRUE, order = c("original", "hclust"), digits = 2,
                           palette = COR_PAL, bin_mark = " [0/1]", title = NULL, subtitle = NULL) {
  order <- match.arg(order)
  o <- .cor_order(cm$r, cm$types, order, split = TRUE)
  R <- cm$r[o, o]; P <- (if (use_adj) cm$p_adj else cm$p)[o, o]; types <- cm$types[o]
  k <- length(o)
  lv <- mark_binary(get_label(o, labels), types, bin_mark)
  fmt <- paste0("%.", digits, "f")
  lab_cols <- paste0("lab__", seq_len(k)); r_cols <- paste0("r__", seq_len(k))
  wide <- data.frame(row_lab = lv, stringsAsFactors = FALSE)
  for (j in seq_len(k)) {
    lab <- ifelse(is.na(R[, j]), "", paste0(sprintf(fmt, R[, j]), p_to_stars(P[, j], ns = "")))
    rr  <- R[, j]
    lab[seq_len(k) < j] <- ""; rr[seq_len(k) < j] <- NA     # upper triangle left empty
    lab[j] <- "1"
    wide[[lab_cols[j]]] <- lab; wide[[r_cols[j]]] <- rr
  }
  meth <- cm$method
  g <- gt::gt(wide, rowname_col = "row_lab")
  if (!is.null(palette)) g <- .gt_colour(g, r_cols, lab_cols, palette)
  g <- g %>%
    gt::cols_hide(r_cols) %>%
    gt::cols_label(.list = stats::setNames(as.list(lv), lab_cols)) %>%
    gt::cols_align("center", columns = lab_cols) %>%
    gt::tab_header(title = title %||% "Correlation matrix",
                   subtitle = subtitle %||% sprintf("%s correlation | pairwise complete N = %s | Computed with %s",
                                                    .meth_name(meth), .n_range_txt(cm$n[upper.tri(cm$n)]),
                                                    .cor_engine_txt(cm$engine)))
  is_bin <- types == "binary"
  if (any(is_bin)) g <- gt::tab_footnote(g, footnote = paste0("Binary variable", bin_mark, ", coded 0/1 (upper level = 1)."),
                                         locations = gt::cells_stub(rows = which(is_bin)))
  g <- .gt_notes(g, list(
    .cor_p_line(if (use_adj) cm$p_adjust else "none", sprintf("the %d distinct pairs", sum(upper.tri(cm$r)))),
    if (order == "hclust") sprintf("Order: hierarchical clustering (Ward D2, distance 1 \u2212 %s), within the continuous and binary blocks.",
                                   .cor_sym_txt(meth)),
    .cor_formula_lines(meth),
    if (any(is_bin) && meth %in% .core_methods) .cor_binary_lines(meth, bin_mark),
    .cor_ref_line()))
  .gt_style(g)
}


#' Table of the effect sizes of binary variables
#'
#' Formats the output of [cor_binary_effects()]: groups and sizes, means (SD),
#' difference in means with its CI and t-test, rank-biserial correlation
#' r_rb with its CI and Mann-Whitney test, and the observed r_s next to its
#' maximum for this split, so that the attenuation of the correlation is visible.
#'
#' @param eff output of [cor_binary_effects()]
#' @param labels named vector variable -> label
#' @param digits decimals
#' @param title,subtitle override the automatic texts
#' @return a `gt_tbl`
#' @export
#' @examples
#' d <- simulate_cor_data(150)
#' tbl_cor_binary(cor_binary_effects(d, c("BIN_CRP", "BIN_IL6"), c("depression", "fatigue")),
#'                labels = attr(d, "labels"))
tbl_cor_binary <- function(eff, labels = NULL, digits = 2, title = NULL, subtitle = NULL) {
  meta <- attr(eff, "cor_meta") %||% list()
  labels <- labels %||% meta$labels
  ci <- meta$ci %||% 0.95
  f <- function(v) ifelse(is.na(v), "", sprintf(paste0("%.", digits, "f"), v))
  fci <- function(lo, hi) ifelse(is.na(lo), "", paste0("[", f(lo), ", ", f(hi), "]"))
  show <- data.frame(
    group  = get_label(eff$Y, labels),
    var    = get_label(eff$X, labels),
    split  = sprintf("%s / %s: %d / %d (%.0f%%)", eff$level0, eff$level1, eff$n0, eff$n1, 100 * eff$prop1),
    m0     = paste0(f(eff$mean0), " (", f(eff$sd0), ")"),
    m1     = paste0(f(eff$mean1), " (", f(eff$sd1), ")"),
    diff   = paste(f(eff$diff), fci(eff$diff_low, eff$diff_high)),
    p_t    = .fmt_p(eff$p_t),
    rb     = paste(f(eff$r_rb), fci(eff$r_rb_low, eff$r_rb_high)),
    rb_num = eff$r_rb,
    p_mw   = .fmt_p(eff$p_mw),
    rs     = paste0(f(eff$r_s), " / ", f(eff$r_s_bound)),
    stringsAsFactors = FALSE)
  lvl <- sprintf("%g%% CI", 100 * ci)
  g <- gt::gt(show, groupname_col = "group") %>%
    .gt_colour("rb_num", "rb") %>%
    gt::cols_hide("rb_num") %>%
    gt::cols_label(var = "Binary variable", split = "Groups 0 / 1: n (% group 1)",
                   m0 = "Mean (SD), group 0", m1 = "Mean (SD), group 1",
                   diff = paste0("M\u2081 \u2212 M\u2080 [", lvl, "]"), p_t = "p (t-test)",
                   rb = paste0("r_rb [", lvl, "]"), p_mw = "p (Mann-Whitney)",
                   rs = "Observed r\u209b / max |r\u209b|") %>%
    gt::tab_spanner("Mean difference", columns = c("diff", "p_t"), id = "sp_diff") %>%
    gt::tab_spanner("Rank-biserial correlation", columns = c("rb", "p_mw"), id = "sp_rb") %>%
    gt::cols_align("center", columns = c("split", "m0", "m1", "diff", "p_t", "rb", "p_mw", "rs")) %>%
    gt::tab_header(title = title %||% "Binary variables: effect sizes to report",
                   subtitle = subtitle %||% sprintf("Group 1 = upper level of each binary variable | N = %s | raw p-values (no adjustment)",
                                                    .n_range_txt(eff$n0 + eff$n1)))
  .gt_style(.gt_notes(g, list(
    sprintf("M\u2081 \u2212 M\u2080: difference in means, Student's t-test with equal variances (= test of the point-biserial r); %s from stats::t.test() {stats}.", lvl),
    paste0("r_rb = 2U / (n\u2080 n\u2081) \u2212 1 = P(Y\u2081 > Y\u2080) \u2212 P(Y\u2081 < Y\u2080) = 2 AUC \u2212 1, U = Mann-Whitney statistic of stats::wilcox.test() {stats} (Kerby 2014, doi:10.2466/11.IT.3.1); ",
           if (identical(meta$rb_ci_engine, "effectsize::rank_biserial")) paste0(lvl, " from effectsize::rank_biserial() {effectsize}.")
           else "CI not computed (package effectsize not installed)."),
    "Observed r\u209b / max |r\u209b|: Spearman correlation with the 0/1 variable and its maximum \u221a(3pq) for this split (no ties): the correlation is bounded by the split and must not be compared with that of a continuous variable (Babchishin & Helmus 2016, doi:10.3758/s13428-015-0627-7).")))
}
