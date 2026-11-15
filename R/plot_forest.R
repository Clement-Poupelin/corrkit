# ====================================================================
# plot_forest.R -- coefficients with their confidence intervals
# ====================================================================

#' Forest plot of a correlation analysis: each coefficient with its CI
#'
#' One row per pair: point = coefficient, bar = confidence interval, filled
#' point = p (adjusted if requested) below `alpha`. A cross analysis has one
#' panel per `y` (e.g. per score); a matrix analysis lists every pair.
#' Complements the heatmap: the heatmap shows the structure, the forest plot
#' shows the precision of each coefficient.
#'
#' @param x output of [cor_compute()]
#' @param labels named vector variable -> label (default: the labels given to [cor_compute()])
#' @param use_adj filled points from adjusted (`TRUE`) or raw p-values
#' @param alpha significance threshold of the filled points
#' @param sort sort the pairs by coefficient within each panel
#' @param base_size base font size
#' @param title title (default "<Method> correlations")
#' @param palette colours of the negative / positive coefficients (named low / high)
#' @param show_note add the binary-variable line to the caption
#' @param caption_width wrap the caption lines at this number of characters (`NULL` = no wrapping)
#' @return a ggplot
#' @seealso `plot()` on a [cor_compute()] result, which calls this function for `type = "forest"`.
#' @export
#' @examples
#' d  <- simulate_cor_data(150)
#' cc <- cor_compute(d, x = c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"),
#'                   y = c("depression", "fatigue", "cognition"), labels = attr(d, "labels"))
#' plot_cor_forest(cc)
plot_cor_forest <- function(x, labels = NULL, use_adj = TRUE, alpha = 0.05, sort = TRUE, base_size = 10,
                            title = NULL, palette = COR_PAL, show_note = TRUE, caption_width = 110) {
  .check_result(x, "plot_cor_forest")
  m <- x$meta; tab <- x$table
  labels <- labels %||% m$labels
  xl <- mark_binary(get_label(tab$X, labels), tab$x_type, m$bin_mark)
  yl <- mark_binary(get_label(tab$Y, labels), tab$y_type, m$bin_mark)
  pv <- if (use_adj) tab$p_adj else tab$p
  d <- data.frame(row = if (m$mode == "matrix") paste(xl, "-", yl) else xl,
                  panel = if (m$mode == "matrix") "" else yl,
                  r = tab$r, lo = tab$ci_low, hi = tab$ci_high,
                  sig = ifelse(!is.na(pv) & pv < alpha, "yes", "no"),
                  dir = ifelse(tab$r >= 0, "pos", "neg"), stringsAsFactors = FALSE)
  d <- d[!is.na(d$r), , drop = FALSE]
  if (m$mode == "cross") d$panel <- factor(d$panel, levels = unique(mark_binary(get_label(m$y_vars, labels),
                                                                                 m$y_types, m$bin_mark)))
  # one key per row and panel, ordered by panel then by coefficient
  d$key <- paste(as.integer(factor(d$panel)), d$row, sep = "@@")
  ord <- if (sort) order(as.integer(factor(d$panel)), d$r) else rev(seq_len(nrow(d)))
  d$key <- factor(d$key, levels = unique(d$key[ord]))
  has_bin <- any(c(tab$x_type, tab$y_type) == "binary")
  sig_lab <- c(yes = sprintf("%s < %g", if (use_adj && !identical(m$p_adjust, "none")) "adjusted p" else "p", alpha),
               no = sprintf("%s >= %g", if (use_adj && !identical(m$p_adjust, "none")) "adjusted p" else "p", alpha))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = r, y = key)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi, colour = dir), width = 0.3,
                           orientation = "y", na.rm = TRUE, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(colour = dir, shape = sig), size = 2.4, fill = "white", stroke = 0.9) +
    ggplot2::scale_colour_manual(values = c(neg = palette[["low"]], pos = palette[["high"]]), guide = "none") +
    ggplot2::scale_shape_manual(values = c(yes = 16, no = 21), labels = sig_lab, name = NULL, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = function(k) sub("^[^@]*@@", "", k)) +
    ggplot2::coord_cartesian(xlim = c(-1, 1)) +
    .labs_safe(title = title %||% sprintf("%s correlations", .meth_name(m$meth)),
               subtitle = .cor_subtitle(m, tab, labels),
               x = sprintf("Coefficient (%s) and %s", .cor_sym_txt(m$meth),
                           if (!is.null(m$ci)) sprintf("%g%% CI", 100 * m$ci) else "no CI"),
               y = NULL,
               caption = .wrap_caption(list(
                 c("How to read: point = coefficient, bar = confidence interval; a bar that does not cross 0 = an association of known sign.",
                   .cor_p_line(if (use_adj) m$p_adjust else "none", if (use_adj) m$family else NULL)),
                 c(.cor_formula_lines(m$meth, m$partial, m$k), .tbl_ci_line(m),
                   if (show_note && has_bin) .cor_binary_lines(m$meth, m$bin_mark) else character(0))),
                 caption_width)) +
    theme_pub(base_size) + .cor_text_theme() +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  if (m$mode == "cross")
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(panel), scales = "free_y", space = "free_y") +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
  p
}
