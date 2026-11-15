# ====================================================================
# plot_heatmap.R -- correlation heatmap of a cor_result (matrix or cross)
# ====================================================================
# Cells involving a binary [0/1] variable are not on the same scale as
# continuous pairs (see cor_binary_note()): they are drawn in a separate block
# (facets or separator lines) and explained in the caption; one colour scale
# for all cells.
# ====================================================================

# Internal: tiles + fill scale; `d` has a column rho
.cor_tiles <- function(d, sym, palette = COR_PAL)
  list(ggplot2::geom_tile(data = d, ggplot2::aes(fill = rho), colour = "white"),
       ggplot2::scale_fill_gradient2(low = palette[["low"]], mid = "white", high = palette[["high"]], midpoint = 0,
                                     limits = c(-1, 1), name = sym))

# Internal: order of the variables of a matrix: "original", or hierarchical
# clustering (Ward D2, distance 1 - r), within the continuous and binary blocks if `split`
.cor_order <- function(R, types, order = "hclust", split = TRUE) {
  vars <- rownames(R)
  is_bin <- types[vars] == "binary"
  split <- split && any(is_bin) && any(!is_bin)
  hc <- function(v) {
    if (order != "hclust" || length(v) <= 2) return(v)
    Rf <- R[v, v, drop = FALSE]; Rf[is.na(Rf)] <- 0
    v[stats::hclust(stats::as.dist(1 - Rf), method = "ward.D2")$order]
  }
  if (split) c(hc(vars[!is_bin]), hc(vars[is_bin])) else hc(vars)
}

#' Heatmap of a correlation analysis
#'
#' * matrix analysis (`cor_compute(y = NULL)`): lower triangle, diagonal included;
#' * cross analysis: one row per `x`, one column per `y`.
#'
#' Each cell shows the coefficient and the stars of its p-value (adjusted if an
#' adjustment was requested); the subtitle gives the method, the adjustment
#' for covariates, N and the engine; the caption gives the status of the
#' p-values, the formula of the coefficient and the reading of binary cells.
#'
#' @param x output of [cor_compute()]
#' @param labels named vector variable -> display label (default: the labels given to [cor_compute()])
#' @param use_adj stars on adjusted (`TRUE`) or raw p-values
#' @param order matrix only: `"hclust"` (within each block if `split_binary`) or `"original"`
#' @param show_values print the coefficient and the stars in each cell
#' @param base_size base font size
#' @param title title (default "Correlation matrix" or "<Method> correlations")
#' @param value_size size of the values
#' @param show_note add the binary-variable line to the caption
#' @param split_binary binary variables in a separate block
#' @param palette colours of the negative / positive correlations (named low / high)
#' @param caption_width wrap the caption lines at this number of characters (`NULL` = no wrapping)
#' @param separator matrix only: `"facet"` (panels with strips, default) or
#'   `"line"` (separator lines, square tiles)
#' @return a ggplot
#' @seealso `plot()` on a [cor_compute()] result, which calls this function for `type = "heatmap"`.
#' @export
#' @examples
#' d  <- simulate_cor_data(150)
#' cm <- cor_compute(d, c("CRP", "IL6", "BDNF", "NfL", "depression", "BIN_CRP", "BIN_IL6"),
#'                   labels = attr(d, "labels"))
#' plot_cor_heatmap(cm)
#' cc <- cor_compute(d, x = c("CRP", "IL6", "NfL", "BIN_CRP"), y = c("depression", "cognition"),
#'                   covariates = "age", labels = attr(d, "labels"))
#' plot_cor_heatmap(cc)
plot_cor_heatmap <- function(x, labels = NULL, use_adj = TRUE, order = c("hclust", "original"),
                             show_values = TRUE, base_size = 10, title = NULL, value_size = NULL,
                             show_note = TRUE, split_binary = TRUE, palette = COR_PAL,
                             caption_width = 110, separator = c("facet", "line")) {
  .check_result(x, "plot_cor_heatmap")
  order <- match.arg(order); separator <- match.arg(separator)
  m <- x$meta
  labels <- labels %||% m$labels
  if (m$mode == "matrix")
    .heatmap_matrix(x, labels, use_adj, order, show_values, base_size, title, value_size %||% 2.6,
                    show_note, split_binary, palette, caption_width, separator)
  else
    .heatmap_cross(x, labels, use_adj, show_values, base_size, title, value_size %||% 2.8,
                   show_note, split_binary, palette, caption_width)
}

# Internal: caption shared by both layouts
.heatmap_caption <- function(m, use_adj, order_txt, has_bin, show_note, caption_width) {
  .wrap_caption(list(
    c(.cor_p_line(if (use_adj) m$p_adjust else "none", if (use_adj) m$family else NULL), order_txt),
    c(.cor_formula_lines(m$meth, m$partial, m$k),
      if (show_note && has_bin) .cor_binary_lines(m$meth, m$bin_mark) else character(0))),
    caption_width)
}

# Internal: lower-triangular heatmap of a matrix analysis
.heatmap_matrix <- function(x, labels, use_adj, order, show_values, base_size, title, value_size,
                            show_note, split_binary, palette, caption_width, separator) {
  m <- x$meta; bin_mark <- m$bin_mark
  R <- as.matrix(x, "r"); P <- as.matrix(x, if (use_adj) "p_adj" else "p")
  types <- m$x_types[rownames(R)]
  is_bin <- types == "binary"
  split  <- split_binary && any(is_bin) && any(!is_bin)
  o <- .cor_order(R, types, order, split_binary)
  R <- R[o, o]; P <- P[o, o]; types <- types[o]
  k  <- length(o); nc <- sum(types != "binary")
  lv <- mark_binary(get_label(o, labels), types, mark = bin_mark)
  sym <- .cor_sym(m$meth)

  long <- tibble::tibble(i = as.vector(row(R)), j = as.vector(col(R)),
                         rho = as.vector(R), p = as.vector(P)) %>%
    dplyr::filter(i >= j) %>%
    dplyr::mutate(x = factor(lv[j], levels = lv), y = factor(lv[i], levels = rev(lv)),
                  lab = dplyr::case_when(is.na(rho) ~ "", i == j ~ "1",
                                         TRUE ~ paste0(sprintf("%.2f", rho), "\n", p_to_stars(p, ns = "", trend = FALSE))))

  order_txt <- if (order == "original") ""
               else sprintf("Order: hierarchical clustering (Ward D2, distance 1 \u2212 %s)%s.",
                            .cor_sym_txt(m$meth), if (split) ", within the continuous and the binary blocks" else "")
  caption <- .heatmap_caption(m, use_adj, order_txt, any(is_bin), show_note, caption_width)

  bin_block <- paste0("Binary", bin_mark)
  long <- long %>%
    dplyr::mutate(x_block = factor(ifelse(types[j] == "binary", bin_block, "Continuous"), levels = c("Continuous", bin_block)),
                  y_block = factor(ifelse(types[i] == "binary", bin_block, "Continuous"), levels = c("Continuous", bin_block)))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = y)) +
    .cor_tiles(long, sym, palette)
  if (show_values) p <- p + ggplot2::geom_text(ggplot2::aes(label = lab), size = value_size, lineheight = 0.8)
  p <- p +
    .labs_safe(title = title %||% "Correlation matrix", subtitle = .cor_subtitle(m, x$table, labels),
               caption = caption, x = NULL, y = NULL) +
    theme_pub(base_size) + rotate_x(45) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) + .cor_text_theme()

  if (split && separator == "facet") {
    # 2 x 2 panels: continuous x continuous (triangle), binary x continuous (rectangle),
    # binary x binary (triangle); the top-right panel (continuous rows x binary columns)
    # is empty and made invisible by removing the panel background and border.
    p + ggplot2::facet_grid(rows = ggplot2::vars(y_block), cols = ggplot2::vars(x_block),
                            scales = "free", space = "free", switch = "x") +
      ggplot2::theme(panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     panel.spacing = ggplot2::unit(0.6, "lines"),
                     strip.placement = "outside")
  } else if (split) {
    p + ggplot2::annotate("segment", x = 0.5, xend = nc + 0.5, y = k - nc + 0.5, yend = k - nc + 0.5,
                          linewidth = 0.9, colour = "grey20") +
      ggplot2::annotate("segment", x = nc + 0.5, xend = nc + 0.5, y = 0.5, yend = k - nc + 0.5,
                        linewidth = 0.9, colour = "grey20") +
      ggplot2::coord_fixed()
  } else {
    p + ggplot2::coord_fixed()
  }
}

# Internal: rectangular heatmap of a cross analysis (rows = x, columns = y)
.heatmap_cross <- function(x, labels, use_adj, show_values, base_size, title, value_size,
                           show_note, split_binary, palette, caption_width) {
  m <- x$meta; tab <- x$table; bin_mark <- m$bin_mark
  x_vars <- m$x_vars; y_vars <- m$y_vars; x_types <- m$x_types; y_types <- m$y_types
  has_bin <- any(c(x_types, y_types) == "binary")
  caption <- .heatmap_caption(m, use_adj, "", has_bin, show_note, caption_width)

  bin_block <- paste0("Binary", bin_mark)
  x_lab <- mark_binary(get_label(x_vars, labels), x_types, mark = bin_mark)
  y_lab <- mark_binary(get_label(y_vars, labels), y_types, mark = bin_mark)
  pv <- if (use_adj) tab$p_adj else tab$p
  d <- tab %>%
    dplyr::mutate(rho = r,
                  Xl = factor(x_lab[match(X, x_vars)], levels = rev(x_lab)),
                  Yl = factor(y_lab[match(Y, y_vars)], levels = y_lab),
                  x_block = factor(ifelse(x_type == "binary", bin_block, "Continuous"), levels = c("Continuous", bin_block)),
                  y_block = factor(ifelse(y_type == "binary", bin_block, "Continuous"), levels = c("Continuous", bin_block)),
                  lab = ifelse(is.na(r), "", paste0(sprintf("%.2f", r), p_to_stars(pv, ns = "", trend = FALSE))))

  pl <- ggplot2::ggplot(d, ggplot2::aes(x = Yl, y = Xl)) +
    .cor_tiles(d, .cor_sym(m$meth), palette)
  if (show_values) pl <- pl + ggplot2::geom_text(ggplot2::aes(label = lab), size = value_size)
  facet_rows <- split_binary && length(unique(x_types == "binary")) == 2
  facet_cols <- split_binary && length(unique(y_types == "binary")) == 2
  if (facet_rows || facet_cols)
    pl <- pl + ggplot2::facet_grid(rows = if (facet_rows) ggplot2::vars(x_block),
                                   cols = if (facet_cols) ggplot2::vars(y_block),
                                   scales = "free", space = "free")
  pl +
    .labs_safe(title = title %||% sprintf("%s correlations", .meth_name(m$meth)),
               subtitle = .cor_subtitle(m, tab, labels), caption = caption, x = NULL, y = NULL) +
    theme_pub(base_size) + rotate_x(45) + ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    .cor_text_theme()
}
