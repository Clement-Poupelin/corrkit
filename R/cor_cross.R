# ====================================================================
# cor_cross.R -- cross correlations X x Y, optionally partial
# ====================================================================

#' Cross correlations X x Y (e.g. biomarkers x scores), optionally partial
#'
#' Every pair (x, y) is tested with [partial_cor_test()]; p-values are then
#' adjusted with [stats::p.adjust()] within the family chosen by `adjust_by`.
#' The result feeds [plot_cor_cross()], [tbl_cor_cross()] and [tbl_cor_long()].
#'
#' @param data data frame
#' @param x_vars predictors (rows of the heatmap, e.g. biomarkers)
#' @param y_vars outcomes (columns, e.g. scores)
#' @param covariates covariates of the partial correlation (`NULL` = unadjusted)
#' @param labels named vector variable -> display label
#' @param p_adjust method of [stats::p.adjust()]
#' @param adjust_by family of the adjustment: `"y"` (the predictors of each
#'   score, default, as in univariable forest plots), `"x"` (the scores of each
#'   predictor), `"all"` (every test), `"none"`
#' @param meth `"spearman"` (default), `"pearson"` or `"kendall"`
#' @param engine `"ppcor"` (default) or `"base"` for the partial correlations
#' @param ci confidence level (`NULL` = none)
#' @param bin_mark marker appended to the labels of binary variables
#' @return object of class `cor_cross`: list(table, meta). `table` has one row
#'   per pair: X, Y, x_type, y_type, rho (the coefficient, whatever the method),
#'   ci_low, ci_high, p, n, k, engine, p_adj.
#' @export
#' @examples
#' d  <- simulate_cor_data(150)
#' cc <- cor_cross(d, x_vars = c("CRP", "IL6", "NfL", "BIN_CRP"),
#'                 y_vars = c("depression", "fatigue", "cognition"),
#'                 covariates = c("age", "sex"), labels = attr(d, "labels"))
#' cc
cor_cross <- function(data, x_vars, y_vars, covariates = NULL, labels = NULL,
                      p_adjust = "BH", adjust_by = c("y", "x", "all", "none"),
                      meth = c("spearman", "pearson", "kendall"), engine = c("ppcor", "base"),
                      ci = 0.95, bin_mark = " [0/1]") {
  meth <- match.arg(meth); adjust_by <- match.arg(adjust_by); engine <- match.arg(engine)
  x_vars <- keep_existing(x_vars, names(data), "cor_cross")
  y_vars <- keep_existing(y_vars, names(data), "cor_cross")
  if (!length(x_vars) || !length(y_vars)) stop("cor_cross: no usable x_vars or y_vars.", call. = FALSE)
  covariates <- if (length(covariates)) keep_existing(covariates, names(data), "cor_cross") else NULL
  if (meth == "kendall" && length(covariates) && !requireNamespace("ppcor", quietly = TRUE))
    stop("cor_cross: partial Kendall correlations need the ppcor package.", call. = FALSE)
  Z <- if (length(covariates)) data[, covariates, drop = FALSE] else NULL
  x_types <- vapply(x_vars, function(v) cor_var_type(data[[v]]), character(1))
  y_types <- vapply(y_vars, function(v) cor_var_type(data[[v]]), character(1))
  all_vars <- unique(c(x_vars, y_vars))
  num <- lapply(stats::setNames(all_vars, all_vars), function(v) as_cor_numeric(data[[v]], v, meth))

  tab <- tidyr::expand_grid(X = x_vars, Y = y_vars) %>%
    dplyr::mutate(res = purrr::map2(X, Y, ~ partial_cor_test(num[[.x]], num[[.y]], Z, meth, engine, ci)),
                  x_type = unname(x_types[X]), y_type = unname(y_types[Y]),
                  rho = purrr::map_dbl(res, "estimate"),
                  ci_low = purrr::map_dbl(res, "ci_low"), ci_high = purrr::map_dbl(res, "ci_high"),
                  p = purrr::map_dbl(res, "p.value"),
                  n = purrr::map_int(res, ~ as.integer(.x$n)),
                  k = purrr::map_int(res, ~ as.integer(.x$k)),
                  engine = purrr::map_chr(res, "engine"),
                  ci_engine = purrr::map_chr(res, "ci_engine")) %>%
    dplyr::select(-res)
  tab <- switch(adjust_by, y = dplyr::group_by(tab, Y), x = dplyr::group_by(tab, X), tab) %>%
    dplyr::mutate(p_adj = if (adjust_by == "none") p else stats::p.adjust(p, p_adjust)) %>%
    dplyr::ungroup()

  family <- switch(adjust_by,
                   y = sprintf("the %d predictors of each score (column)", length(x_vars)),
                   x = sprintf("the %d scores of each predictor (row)", length(y_vars)),
                   all = sprintf("all %d tests", nrow(tab)),
                   none = NULL)
  partial <- length(covariates) > 0
  meta <- list(meth = meth, p_adjust = if (adjust_by == "none") "none" else p_adjust, adjust_by = adjust_by,
               family = family, covariates = covariates, partial = partial,
               k = if (partial) suppressWarnings(max(tab$k, na.rm = TRUE)) else 0L,
               x_vars = x_vars, y_vars = y_vars, x_types = x_types, y_types = y_types,
               labels = labels, bin_mark = bin_mark, ci = ci,
               engine = unique(tab$engine), ci_engine = setdiff(unique(tab$ci_engine), "none"))
  attr(tab, "cor_meta") <- meta
  structure(list(table = tab, meta = meta), class = c("cor_cross", "list"))
}

# Internal: accept a cor_cross object, the list returned by plot_cor_cross(), or its $table
.as_cor_cross <- function(x) {
  if (inherits(x, "cor_cross")) return(x)
  tab <- if (is.data.frame(x)) x else if (is.list(x) && is.data.frame(x$table)) x$table else NULL
  meta <- if (!is.null(tab)) attr(tab, "cor_meta") else NULL
  if (is.null(meta)) stop("Expected the output of cor_cross() or plot_cor_cross().", call. = FALSE)
  structure(list(table = tab, meta = meta), class = c("cor_cross", "list"))
}

# Internal: subtitle shared by the plot and the tables of a cor_cross
.cor_cross_subtitle <- function(meta, tab, labels = meta$labels, newline = "\n") {
  sprintf("%s %s | N = %s%sComputed with %s; p adjustment: %s",
          .meth_name(meta$meth),
          if (meta$partial) paste0("partial, adjusted for ", paste(get_label(meta$covariates, labels), collapse = ", "))
          else "unadjusted",
          .n_range_txt(tab$n), newline, .cor_engine_txt(tab$engine),
          if (identical(meta$p_adjust, "none")) "none" else "stats::p.adjust() {stats}")
}

#' @export
print.cor_cross <- function(x, digits = 2, ...) {
  m <- x$meta; tab <- x$table
  cat(sprintf("<cor_cross> %d x %d %s correlations (%s)\n", length(m$x_vars), length(m$y_vars), .meth_name(m$meth),
              if (m$partial) paste("partial, adjusted for", paste(m$covariates, collapse = ", ")) else "unadjusted"))
  cat("  ", .cor_p_line(m$p_adjust, m$family), "\n", sep = "")
  wide <- tab %>%
    dplyr::mutate(lab = ifelse(is.na(rho), "", paste0(sprintf(paste0("%.", digits, "f"), rho),
                                                     p_to_stars(p_adj, ns = "")))) %>%
    dplyr::select(X, Y, lab) %>%
    tidyr::pivot_wider(names_from = Y, values_from = lab)
  print(as.data.frame(wide), row.names = FALSE)
  invisible(x)
}
