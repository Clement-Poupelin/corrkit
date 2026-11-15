# ====================================================================
# methods.R -- plot(), autoplot() and summary() of a cor_result
# ====================================================================
# Conventions followed (as in correlation / see, corrr, broom):
#   print()   short display of the result (in cor_compute.R)
#   summary() characteristics of the analysis and content of the object
#   plot()    a ggplot, whose kind is chosen with `type`
#   autoplot() the same ggplot (ggplot2 generic)
#   as.matrix(), as.data.frame() the numbers, for further work
#   tbl_cor() the gt tables (a function rather than a generic: dplyr::tbl()
#             already exists with another meaning)
# ====================================================================

#' Plot a correlation analysis
#'
#' `plot()` (and `autoplot()`, its ggplot2 equivalent) draw the output of
#' [cor_compute()]; the kind of plot is chosen with `type`:
#' * `"heatmap"` (default): [plot_cor_heatmap()], triangle for a matrix,
#'   rectangle for a cross analysis;
#' * `"forest"`: [plot_cor_forest()], each coefficient with its confidence interval;
#' * `"scatter"`: [plot_cor_pair()] for one pair (`pair = c("x", "y")`; by
#'   default the pair with the largest |coefficient|), with the covariates and
#'   the method of the analysis.
#'
#' @param x,object output of [cor_compute()]
#' @param type `"heatmap"`, `"forest"` or `"scatter"`
#' @param pair `type = "scatter"`: the two variables of the pair
#' @param ... passed to the plotting function (`labels`, `title`, `base_size`, `palette`...)
#' @return a ggplot
#' @export
#' @examples
#' d   <- simulate_cor_data(150)
#' res <- cor_compute(d, x = c("CRP", "IL6", "NfL", "BIN_CRP"), y = c("depression", "fatigue"),
#'                    covariates = "age", labels = attr(d, "labels"))
#' plot(res)
#' plot(res, type = "forest")
#' plot(res, type = "scatter", pair = c("IL6", "fatigue"))
plot.cor_result <- function(x, type = c("heatmap", "forest", "scatter"), pair = NULL, ...) {
  type <- match.arg(type)
  switch(type,
         heatmap = plot_cor_heatmap(x, ...),
         forest  = plot_cor_forest(x, ...),
         scatter = .plot_scatter(x, pair, ...))
}

#' @rdname plot.cor_result
#' @importFrom ggplot2 autoplot
#' @export
autoplot.cor_result <- function(object, type = c("heatmap", "forest", "scatter"), pair = NULL, ...)
  plot.cor_result(object, type = type, pair = pair, ...)

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

# Internal: scatter plot of one pair of a cor_result
.plot_scatter <- function(x, pair = NULL, ...) {
  m <- x$meta; tab <- x$table
  if (is.null(pair)) {
    best <- which.max(abs(tab$r))
    if (!length(best)) stop("plot: no pair with a coefficient.", call. = FALSE)
    pair <- c(tab$X[best], tab$Y[best])
    message(sprintf("plot(type = \"scatter\"): pair with the largest |coefficient|: %s - %s", pair[1], pair[2]))
  }
  if (length(pair) != 2 || !all(pair %in% names(x$data)))
    stop("plot: `pair` must name two variables of the analysis.", call. = FALSE)
  meth <- if (m$meth %in% .core_methods) m$meth else {
    message(sprintf("plot(type = \"scatter\"): method \"%s\" not available for one pair; using spearman.", m$meth))
    "spearman"
  }
  args <- list(...)
  defaults <- list(data = x$data, x = pair[1], y = pair[2], covariates = m$covariates,
                   labels = m$labels, meth = meth, ci = m$ci %||% 0.95, bin_mark = m$bin_mark,
                   engine = if ("residuals (base R)" %in% m$engine) "base" else "ppcor")
  do.call(plot_cor_pair, utils::modifyList(defaults, args))
}


#' Summary of a correlation analysis
#'
#' Characteristics of the analysis (mode, method, covariates, engines,
#' confidence intervals, adjustment of the p-values), the variables and their
#' types, the number of pairs and of significant pairs, the strongest pairs,
#' and the content of the object.
#'
#' @param object output of [cor_compute()]
#' @param top number of strongest pairs to list
#' @param alpha significance threshold used in the counts
#' @param ... unused
#' @return object of class `summary.cor_result`: list(settings, variables,
#'   pairs, strongest, contents), printed by `print()`
#' @export
#' @examples
#' d <- simulate_cor_data(150)
#' summary(cor_compute(d, c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"), labels = attr(d, "labels")))
summary.cor_result <- function(object, top = 5, alpha = 0.05, ...) {
  m <- object$meta; tab <- object$table
  all_vars <- unique(c(m$x_vars, m$y_vars))
  role <- ifelse(all_vars %in% m$x_vars & all_vars %in% m$y_vars & m$mode == "cross", "x and y",
                 ifelse(all_vars %in% m$x_vars, if (m$mode == "matrix") "matrix" else "x", "y"))
  types <- c(m$x_types, m$y_types)[all_vars]
  variables <- data.frame(variable = all_vars, label = get_label(all_vars, m$labels), role = role,
                          type = unname(types),
                          n = vapply(all_vars, function(v) sum(!is.na(object$data[[v]])), numeric(1)),
                          row.names = NULL, stringsAsFactors = FALSE)
  if (length(m$covariates))
    variables <- rbind(variables, data.frame(
      variable = m$covariates, label = get_label(m$covariates, m$labels), role = "covariate",
      type = vapply(m$covariates, function(v) cor_var_type(object$data[[v]]), character(1)),
      n = vapply(m$covariates, function(v) sum(!is.na(object$data[[v]])), numeric(1)),
      row.names = NULL, stringsAsFactors = FALSE))
  adjusted <- !identical(m$p_adjust, "none")
  pairs <- list(n_pairs = nrow(tab), n_range = .n_range_txt(tab$n),
                n_binary = sum(tab$x_type == "binary" | tab$y_type == "binary"),
                n_na = sum(is.na(tab$r)),
                sig_raw = sum(tab$p < alpha, na.rm = TRUE),
                sig_adj = if (adjusted) sum(tab$p_adj < alpha, na.rm = TRUE) else NA_integer_,
                alpha = alpha)
  o <- order(-abs(tab$r))
  strongest <- as.data.frame(tab[utils::head(o[!is.na(tab$r[o])], top),
                                 c("X", "Y", "r", "ci_low", "ci_high", "p", "p_adj", "n")])
  settings <- list(
    mode = if (m$mode == "matrix") sprintf("matrix: every pair of %d variables", length(m$x_vars))
           else sprintf("cross: %d x %d variables", length(m$x_vars), length(m$y_vars)),
    method = .meth_name(m$meth),
    covariates = if (m$partial) sprintf("%s (%d term%s)", paste(m$covariates, collapse = ", "), m$k,
                                        if (m$k > 1) "s" else "") else "none (unadjusted)",
    engine = .cor_engine_txt(tab$engine),
    ci = if (is.null(m$ci)) "none" else sprintf("%g%% (%s)", 100 * m$ci,
                                               if (length(m$ci_engine)) paste(m$ci_engine, collapse = ", ") else "none"),
    p_values = if (adjusted) sprintf("adjusted, %s over %s", .p_adjust_name(m$p_adjust), m$family)
               else "raw (no multiple-testing adjustment)")
  contents <- c(
    "$table" = sprintf("one row per pair (%d x %d): %s", nrow(tab), ncol(tab), paste(names(tab), collapse = ", ")),
    "$meta"  = sprintf("settings of the analysis (%s)", paste(names(m), collapse = ", ")),
    "$data"  = sprintf("the %d columns used (%d rows)", ncol(object$data), nrow(object$data)))
  structure(list(settings = settings, variables = variables, pairs = pairs, strongest = strongest,
                 contents = contents, labels = m$labels, sym = .cor_sym_txt(m$meth), adjusted = adjusted),
            class = "summary.cor_result")
}

#' @export
print.summary.cor_result <- function(x, digits = 2, ...) {
  s <- x$settings; p <- x$pairs
  cat("== Correlation analysis (cor_result) ==\n\n")
  cat(sprintf("  %-12s %s\n", c("Mode", "Method", "Covariates", "Engine", "CI", "p-values"),
              c(s$mode, s$method, s$covariates, s$engine, s$ci, s$p_values)), sep = "")
  cat("\n-- Variables --\n")
  print(x$variables, row.names = FALSE)
  cat(sprintf("\n-- Pairs --\n  %d pairs (%d involving a binary variable), N = %s%s\n",
              p$n_pairs, p$n_binary, p$n_range,
              if (p$n_na) sprintf(", %d without a coefficient", p$n_na) else ""))
  cat(sprintf("  p < %g: %d (raw)%s\n", p$alpha, p$sig_raw,
              if (x$adjusted) sprintf(", %d (adjusted)", p$sig_adj) else ""))
  if (nrow(x$strongest)) {
    cat(sprintf("\n-- Strongest pairs (|%s|) --\n", x$sym))
    st <- x$strongest
    out <- data.frame(pair = paste(get_label(st$X, x$labels), "-", get_label(st$Y, x$labels)),
                      r = round(st$r, digits),
                      CI = ifelse(is.na(st$ci_low), "",
                                  sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), st$ci_low, st$ci_high)),
                      p = .fmt_p(st$p), n = st$n, stringsAsFactors = FALSE)
    if (x$adjusted) out$p_adj <- .fmt_p(st$p_adj)
    names(out)[2] <- x$sym
    print(out, row.names = FALSE)
  }
  cat("\n-- Content of the object --\n")
  cat(sprintf("  %-7s %s\n", names(x$contents), x$contents), sep = "")
  cat("  Use: plot(x, type = \"heatmap\" | \"forest\" | \"scatter\"), tbl_cor(x, layout = ...),\n",
      "       as.matrix(x, what = \"r\" | \"p\" | \"p_adj\" | \"n\" | \"ci_low\" | \"ci_high\"), as.data.frame(x)\n", sep = "")
  invisible(x)
}
