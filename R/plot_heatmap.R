# ====================================================================
# plot_heatmap.R -- correlation heatmaps (matrix and cross)
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

#' Lower-triangular heatmap of a correlation matrix
#'
#' @param cm output of [cor_test_matrix()]
#' @param labels named vector variable -> display label
#' @param use_adj stars on adjusted (`TRUE`) or raw p-values
#' @param order `"hclust"` (within each block if `split_binary`) or `"original"`
#' @param show_values print the coefficient and the stars in each cell
#' @param base_size base font size
#' @param title title (default "Correlation matrix")
#' @param value_size size of the values
#' @param show_note add the binary-variable line to the caption
#' @param bin_mark marker appended to the labels of binary variables
#' @param split_binary binary variables in a separate block (last rows / columns)
#' @param palette colours of the negative / positive correlations (named low / high)
#' @param caption_width wrap the caption lines at this number of characters (`NULL` = no wrapping)
#' @param separator `"facet"` (2 x 2 panels with strips, default) or `"line"` (separator lines, square tiles)
#' @return a ggplot
#' @export
#' @examples
#' d  <- simulate_cor_data(150)
#' cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "NfL", "depression", "BIN_CRP", "BIN_IL6"))
#' plot_cor_heatmap(cm, labels = attr(d, "labels"))
plot_cor_heatmap <- function(cm, labels = NULL, use_adj = TRUE, order = c("hclust", "original"),
                             show_values = TRUE, base_size = 10, title = NULL, value_size = 2.6,
                             show_note = TRUE, bin_mark = " [0/1]", split_binary = TRUE,
                             palette = COR_PAL, caption_width = 110,
                             separator = c("facet", "line")) {
  order <- match.arg(order); separator <- match.arg(separator)
  R <- cm$r; P <- if (use_adj) cm$p_adj else cm$p
  vars  <- rownames(R)
  types <- (cm$types %||% stats::setNames(rep("continuous", length(vars)), vars))[vars]
  is_bin <- types == "binary"
  split  <- split_binary && any(is_bin) && any(!is_bin)
  o <- .cor_order(R, types, order, split_binary)
  R <- R[o, o]; P <- P[o, o]; types <- types[o]
  k  <- length(o); nc <- sum(types != "binary")
  lv <- mark_binary(get_label(o, labels), types, mark = bin_mark)
  meth <- cm$method %||% "spearman"
  sym  <- .cor_sym(meth)

  long <- tibble::tibble(i = as.vector(row(R)), j = as.vector(col(R)),
                         rho = as.vector(R), p = as.vector(P)) %>%
    dplyr::filter(i >= j) %>%
    dplyr::mutate(x = factor(lv[j], levels = lv), y = factor(lv[i], levels = rev(lv)),
                  lab = dplyr::case_when(is.na(rho) ~ "", i == j ~ "1",
                                         TRUE ~ paste0(sprintf("%.2f", rho), "\n", p_to_stars(p, ns = "", trend = FALSE))))

  order_txt <- if (order == "original") ""
               else sprintf("Order: hierarchical clustering (Ward D2, distance 1 \u2212 %s)%s.",
                            .cor_sym_txt(meth), if (split) ", within the continuous and the binary blocks" else "")
  n_pairs <- sum(upper.tri(cm$r))
  caption <- .wrap_caption(list(
    c(.cor_p_line(if (use_adj) cm$p_adjust else "none", sprintf("the %d distinct pairs", n_pairs)), order_txt),
    c(.cor_formula_lines(meth),
      if (show_note && any(is_bin)) .cor_binary_lines(meth, bin_mark) else character(0))),
    caption_width)

  bin_block <- paste0("Binary", bin_mark)
  long <- long %>%
    dplyr::mutate(x_block = factor(ifelse(types[j] == "binary", bin_block, "Continuous"), levels = c("Continuous", bin_block)),
                  y_block = factor(ifelse(types[i] == "binary", bin_block, "Continuous"), levels = c("Continuous", bin_block)))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = y)) +
    .cor_tiles(long, sym, palette)
  if (show_values) p <- p + ggplot2::geom_text(ggplot2::aes(label = lab), size = value_size, lineheight = 0.8)
  p <- p +
    ggplot2::labs(title = title %||% "Correlation matrix",
                  subtitle = sprintf("%s correlation | pairwise complete N = %s\nComputed with %s; p adjustment: %s",
                                     .meth_name(meth), .n_range_txt(cm$n[upper.tri(cm$n)]),
                                     .cor_engine_txt(cm$engine %||% "stats::cor.test"),
                                     if (use_adj && !identical(cm$p_adjust, "none")) "stats::p.adjust() {stats}" else "none"),
                  caption = caption, x = NULL, y = NULL) +
    theme_pub(base_size) + rotate_x(45) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) + .cor_text_theme()

  if (split && separator == "facet") {
    # 2 x 2 panels: continuous x continuous (triangle), binary x continuous (rectangle),
    # binary x binary (triangle); the top-right panel (continuous rows x binary columns)
    # is empty and made invisible by removing the panel background and border.
    # coord_fixed() is not allowed with free scales: square tiles come from space = "free" + a square figure.
    p + ggplot2::facet_grid(rows = ggplot2::vars(y_block), cols = ggplot2::vars(x_block),
                            scales = "free", space = "free", switch = "x") +
      ggplot2::theme(panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     panel.spacing = ggplot2::unit(0.6, "lines"),
                     strip.placement = "outside")
  } else if (split) {
    # thin separator lines between the continuous block (top-left) and the binary block
    p + ggplot2::annotate("segment", x = 0.5, xend = nc + 0.5, y = k - nc + 0.5, yend = k - nc + 0.5,
                          linewidth = 0.9, colour = "grey20") +
      ggplot2::annotate("segment", x = nc + 0.5, xend = nc + 0.5, y = 0.5, yend = k - nc + 0.5,
                        linewidth = 0.9, colour = "grey20") +
      ggplot2::coord_fixed()
  } else {
    p + ggplot2::coord_fixed()
  }
}


#' Heatmap of cross correlations X x Y, optionally partial
#'
#' Either computes the correlations ([cor_cross()]) or draws an existing
#' `cor_cross` object passed as `data`.
#'
#' @inheritParams cor_cross
#' @param data data frame, or the output of [cor_cross()]
#' @param base_size base font size
#' @param title title (default "<Method> correlations")
#' @param value_size size of the values
#' @param show_note add the short reading note for binary variables to the caption
#' @param split_binary binary variables in a separate facet ("Continuous" / "Binary [0/1]")
#' @param palette colours of the negative / positive correlations (named low / high)
#' @param caption_width wrap the caption lines at this number of characters
#'   (`NULL` = no wrapping); ~60 for two heatmaps side by side
#' @return list(table, plot, note): `table` = the table of [cor_cross()]
#'   (usable by [tbl_cor_cross()] and [tbl_cor_long()]); `note` = full reading
#'   note for binary variables (`NULL` if none)
#' @export
#' @examples
#' d <- simulate_cor_data(150)
#' res <- plot_cor_cross(d, x_vars = c("CRP", "IL6", "NfL", "BIN_CRP", "BIN_IL6"),
#'                       y_vars = c("depression", "fatigue", "cognition"),
#'                       covariates = c("age", "sex"), labels = attr(d, "labels"))
#' res$plot
plot_cor_cross <- function(data, x_vars, y_vars, covariates = NULL, labels = NULL,
                           p_adjust = "BH", adjust_by = c("y", "x", "all", "none"),
                           meth = c("spearman", "pearson", "kendall"), engine = c("ppcor", "base"),
                           base_size = 10, title = NULL, value_size = 2.8, show_note = TRUE,
                           bin_mark = " [0/1]", split_binary = TRUE,
                           palette = COR_PAL, caption_width = 90, ci = 0.95) {
  cc <- if (inherits(data, "cor_cross")) data
        else cor_cross(data, x_vars, y_vars, covariates, labels, p_adjust, adjust_by, meth, engine, ci, bin_mark)
  m <- cc$meta; tab <- cc$table
  labels <- labels %||% m$labels
  meth <- m$meth
  x_vars <- m$x_vars; y_vars <- m$y_vars; x_types <- m$x_types; y_types <- m$y_types
  sym <- .cor_sym(meth)
  has_bin <- any(c(x_types, y_types) == "binary")

  caption <- .wrap_caption(list(
    .cor_p_line(m$p_adjust, m$family),
    c(.cor_formula_lines(meth, m$partial, m$k),
      if (show_note && has_bin) .cor_binary_lines(meth, bin_mark) else character(0))),
    caption_width)

  bin_block <- paste0("Binary", bin_mark)
  x_lab <- mark_binary(get_label(x_vars, labels), x_types, mark = bin_mark)
  y_lab <- mark_binary(get_label(y_vars, labels), y_types, mark = bin_mark)
  d <- tab %>%
    dplyr::mutate(Xl = factor(x_lab[match(X, x_vars)], levels = rev(x_lab)),
                  Yl = factor(y_lab[match(Y, y_vars)], levels = y_lab),
                  x_block = factor(ifelse(x_type == "binary", bin_block, "Continuous"), levels = c("Continuous", bin_block)),
                  y_block = factor(ifelse(y_type == "binary", bin_block, "Continuous"), levels = c("Continuous", bin_block)),
                  lab = ifelse(is.na(rho), "", paste0(sprintf("%.2f", rho), p_to_stars(p_adj, ns = "", trend = FALSE))))

  pl <- ggplot2::ggplot(d, ggplot2::aes(x = Yl, y = Xl)) +
    .cor_tiles(d, sym, palette) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = value_size)
  facet_rows <- split_binary && length(unique(x_types == "binary")) == 2
  facet_cols <- split_binary && length(unique(y_types == "binary")) == 2
  if (facet_rows || facet_cols)
    pl <- pl + ggplot2::facet_grid(rows = if (facet_rows) ggplot2::vars(x_block),
                                   cols = if (facet_cols) ggplot2::vars(y_block),
                                   scales = "free", space = "free")
  pl <- pl +
    ggplot2::labs(title = title %||% sprintf("%s correlations", .meth_name(meth)),
                  subtitle = .cor_cross_subtitle(m, tab, labels),
                  caption = caption, x = NULL, y = NULL) +
    theme_pub(base_size) + rotate_x(45) + ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    .cor_text_theme()
  list(table = tab, plot = pl, note = if (has_bin) cor_binary_note(meth, "en", mark = bin_mark) else NULL)
}
