# ====================================================================
# utils.R -- small helpers used by every other file
# ====================================================================

#' Default value when NULL
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Keep the variables present in the data, with a warning for the others
#' @noRd
keep_existing <- function(vars, available, where = "corrkit") {
  vars <- unique(as.character(vars))
  miss <- setdiff(vars, available)
  if (length(miss))
    warning(sprintf("%s: variable(s) not found and ignored: %s", where, paste(miss, collapse = ", ")),
            call. = FALSE)
  intersect(vars, available)
}

#' Display label of each variable: `labels[var]` if given, else the variable name
#' @param labels named character vector (or named list) variable -> label
#' @noRd
get_label <- function(vars, labels = NULL) {
  vars <- as.character(vars)
  if (is.null(labels) || !length(vars)) return(vars)
  labels <- unlist(labels)
  out <- unname(labels[vars])
  out[is.na(out)] <- vars[is.na(out)]
  as.character(out)
}

#' Significance stars of p-values
#'
#' `***` p < 0.001, `**` p < 0.01, `*` p < 0.05; otherwise `ns`
#' (or `.` for p < 0.10 when `trend = TRUE`). Missing p -> "".
#' @param p numeric vector of p-values
#' @param ns text for a non-significant p
#' @param trend mark 0.05 <= p < 0.10 with "."
#' @return character vector, same length as `p`
#' @export
#' @examples
#' p_to_stars(c(0.0004, 0.003, 0.04, 0.07, 0.5, NA), ns = "", trend = TRUE)
p_to_stars <- function(p, ns = "ns", trend = FALSE) {
  out <- rep(ns, length(p))
  if (trend) out[!is.na(p) & p < 0.10] <- "."
  out[!is.na(p) & p < 0.05]  <- "*"
  out[!is.na(p) & p < 0.01]  <- "**"
  out[!is.na(p) & p < 0.001] <- "***"
  out[is.na(p)] <- ""
  out
}

#' ggplot2 theme of the corrkit figures
#'
#' Minimal theme; title and caption aligned on the whole plot, caption (the
#' reading guide) left-aligned in grey.
#' @param base_size base font size
#' @export
theme_pub <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.title.position = "plot", plot.caption.position = "plot",
                   plot.subtitle = ggplot2::element_text(colour = "grey25", size = ggplot2::rel(0.9)),
                   plot.caption = ggplot2::element_text(hjust = 0, colour = "grey35", size = ggplot2::rel(0.8)),
                   strip.text = ggplot2::element_text(face = "bold"))
}

#' Rotated x-axis labels
#' @noRd
rotate_x <- function(angle = 45)
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = 1, vjust = 1))

#' p-value as text: "< 0.001" or 3 decimals
#' @noRd
.fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))

#' "p = 0.012" / "p < 0.001"
#' @noRd
.fmt_p_eq <- function(p) ifelse(is.na(p), "p = NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))

#' Range of sample sizes as text: "98" or "95-100"
#' @noRd
.n_range_txt <- function(n) {
  n <- n[is.finite(n)]
  if (!length(n)) return("NA")
  if (min(n) == max(n)) as.character(min(n)) else paste0(min(n), "\u2013", max(n))
}

#' Options of the plot texts
#'
#' * `corrkit.unicode` (default `FALSE`): keep the mathematical characters
#'   (square root, minus sign, tau, less-or-equal...) in the plot texts. The
#'   pdf() and postscript() devices, used by R CMD check and by
#'   `ggsave("x.pdf")`, cannot draw them; set it to `TRUE` if you save your
#'   figures as PNG or with `ggsave(..., device = cairo_pdf)`.
#' * `corrkit.markdown` (default `TRUE`): when the ggtext package is
#'   installed, subtitles, captions and axis titles are rendered as markdown,
#'   so that the notation is typeset: r_s becomes an italic r with a subscript s.
#' @noRd
.use_unicode  <- function() isTRUE(getOption("corrkit.unicode", FALSE))
.use_markdown <- function() isTRUE(getOption("corrkit.markdown", TRUE)) && requireNamespace("ggtext", quietly = TRUE)

#' Text of a ggplot restricted to Latin-1 characters (unless option corrkit.unicode = TRUE)
#'
#' The pdf() and postscript() devices (used by R CMD check, and by ggsave(*.pdf))
#' cannot draw characters outside Latin-1 such as the square root or the
#' minus sign; those are replaced by plain equivalents. Tables (HTML) keep them.
#' @noRd
.plot_text <- function(x) {
  if (is.null(x) || .use_unicode()) return(x)
  map <- c("\u221a" = "sqrt", "\u2212" = "-", "\u2013" = "-", "\u2014" = "-", "\u2261" = "=",
           "\u2264" = "<=", "\u2265" = ">=", "\u2248" = "~", "\u03c4" = "tau", "\u03c6" = "phi",
           "\u03c1" = "rho", "\u2080" = "0", "\u2081" = "1", "\u209b" = "s", "\u2192" = "->")
  for (k in names(map)) x <- gsub(k, map[[k]], x, fixed = TRUE)
  x
}

#' Plain text with the corrkit notation -> markdown/HTML for ggtext
#'
#' r_s, r_xy.z, e_x.z, r_rb, M_1, SD_y, n_0, rank_x, tau_b... become an italic
#' symbol with a subscript (as LaTeX would typeset r_{s}); only these bases are
#' converted, so that function names such as cor_to_ci() are left as they are.
#' Then <, >, & are escaped as HTML entities, the markdown characters * _ [ ]
#' are escaped, and line breaks become <br>.
#' @noRd
.md_math <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  sub_pat <- function(bases) paste0("(?<![A-Za-z0-9_:.])(", paste(bases, collapse = "|"),
                                    ")_([A-Za-z0-9]+(?:[.,][A-Za-z0-9]+)*)")
  # operators upright (rank_x, SD_y), symbols in italics (r_s, e_x.z, M_1, n_0, tau_b)
  x <- gsub(sub_pat(c("rank", "SD")), "\\1<sub>\\2</sub>", x, perl = TRUE)
  x <- gsub(sub_pat(c("tau", "\u03c4", "r", "e", "M", "n")), "<i>\\1</i><sub>\\2</sub>", x, perl = TRUE)
  x <- gsub("([*_\\[\\]])", "\\\\\\1", x, perl = TRUE)
  gsub("\n", "<br>", x, fixed = TRUE)
}

#' The two levels of a binary variable, lower level first (coded 0, upper level coded 1)
#' @noRd
.binary_levels <- function(x) {
  if (is.logical(x)) return(c("FALSE", "TRUE"))
  if (is.factor(x) || is.character(x)) return(levels(droplevels(as.factor(x))))
  as.character(sort(unique(stats::na.omit(x))))
}

#' ggplot2::labs() for the corrkit plots: title, subtitle, caption and axis
#' titles made Latin-1-safe (.plot_text()) and, with ggtext, typeset (.md_math());
#' the matching theme elements come from .cor_text_theme()
#' @noRd
.labs_safe <- function(...) {
  a <- list(...)
  md <- .use_markdown()
  for (k in intersect(c("title", "subtitle", "caption", "x", "y"), names(a))) {
    if (!is.character(a[[k]])) next
    a[[k]] <- .plot_text(a[[k]])
    if (md && k != "title") a[[k]] <- .md_math(a[[k]])
  }
  do.call(ggplot2::labs, a)
}
