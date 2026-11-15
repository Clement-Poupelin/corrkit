# ====================================================================
# cor_compute.R -- the single computing function of corrkit
# ====================================================================
# One operation behind every correlation analysis: a set of pairs (x, y),
# possibly adjusted for covariates.
#   y = NULL : every pair of `x` (matrix, i < j)
#   y given  : every pair x x y (cross)
# The engines are existing functions (stats, ppcor, psych, correlation);
# p-values are raw by default and can be adjusted over an explicit family.
# The result (class "cor_result") feeds plot(), autoplot(), tbl_cor(),
# as.matrix(), as.data.frame() and cor_binary_effects().
# ====================================================================

#' Compute correlations: matrix, cross or single pair, simple or partial
#'
#' `cor_compute()` is the single computing function of corrkit.
#' * `y = NULL`: every pair of the variables `x` (correlation matrix);
#' * `y` given: every pair `x` x `y` (e.g. biomarkers x scores);
#' * one `x` and one `y`: a single pair;
#' * `covariates`: partial correlations (any of the above).
#'
#' **Engines** (an existing function for every computation)
#' * matrix without covariates: `psych::corr.test(adjust = "none")` (Pearson,
#'   Spearman), `correlation::correlation(p_adjust = "none")` (Kendall and
#'   the robust methods), or [stats::cor.test()] pair by pair;
#' * other cases: [partial_cor_test()], i.e. [stats::cor.test()] without
#'   covariates and `ppcor::pcor.test()` with covariates; robust methods use
#'   `correlation::cor_test()`.
#'
#' **p-values** are raw by default (`p_adjust = "none"`). Otherwise they are
#' adjusted with [stats::p.adjust()] over the k(k - 1)/2 pairs (matrix) or
#' over the family set by `adjust_by` (cross).
#'
#' **Binary variables** are coded 0/1, the upper level being 1 (see
#' [as_cor_numeric()]), and flagged in the result; see [cor_binary_note()]
#' and [cor_binary_effects()].
#'
#' @param data data frame
#' @param x variables (rows of the outputs)
#' @param y variables crossed with `x` (columns); `NULL` = every pair of `x`
#' @param covariates covariates of a partial correlation (`NULL` = none);
#'   factors allowed (dummy coded)
#' @param method `"spearman"` (default), `"pearson"`, `"kendall"`; without
#'   covariates also the methods of the correlation package: `"biweight"`,
#'   `"percentage"`, `"blomqvist"`, `"hoeffding"`, `"gaussian"`, `"shepherd"`,
#'   `"distance"`
#' @param p_adjust method of [stats::p.adjust()]: `"none"` (default, raw
#'   p-values) or e.g. `"BH"`, `"holm"`, `"bonferroni"`
#' @param adjust_by family of the adjustment of a cross analysis: `"y"` (the
#'   predictors of each column, default), `"x"` (the columns of each
#'   predictor) or `"all"`; a matrix is always adjusted over all its pairs
#' @param ci confidence level of the intervals (`NULL` = none)
#' @param engine `"auto"` (default), `"psych"`, `"correlation"`, `"stats"`
#'   (matrix without covariates), `"ppcor"` or `"base"` (partial correlations)
#' @param labels named vector variable -> display label, used by the outputs
#' @param bin_mark marker appended to the labels of binary variables
#' @return object of class `cor_result`: a list with
#'   * `table`: one row per pair, with X, Y, x_type, y_type, r (the
#'     coefficient, whatever the method), ci_low, ci_high, p, p_adj, n, k
#'     (number of covariate terms), engine;
#'   * `meta`: settings of the analysis (mode, method, adjustment, family,
#'     covariates, variables and their types, labels...);
#'   * `data`: the columns used, for [plot_cor_pair()] and [cor_binary_effects()].
#'
#'   Methods: `print()`, `plot()`, `autoplot()`, [tbl_cor()], `as.matrix()`,
#'   `as.data.frame()`.
#' @export
#' @examples
#' d   <- simulate_cor_data(150)
#' lab <- attr(d, "labels")
#'
#' # every pair of a set of variables
#' cm <- cor_compute(d, x = c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"), labels = lab)
#' cm
#' as.matrix(cm)
#'
#' # biomarkers x scores, adjusted for age and sex, BH over the predictors of each score
#' cc <- cor_compute(d, x = c("CRP", "IL6", "NfL", "BIN_CRP"), y = c("depression", "cognition"),
#'                   covariates = c("age", "sex"), p_adjust = "BH", labels = lab)
#' cc
cor_compute <- function(data, x, y = NULL, covariates = NULL, method = "spearman",
                        p_adjust = "none", adjust_by = c("y", "x", "all"), ci = 0.95,
                        engine = c("auto", "psych", "correlation", "stats", "ppcor", "base"),
                        labels = NULL, bin_mark = " [0/1]") {
  method    <- match.arg(method, c(.core_methods, .extra_methods))
  adjust_by <- match.arg(adjust_by)
  engine    <- match.arg(engine)
  if (!is.data.frame(data)) stop("cor_compute: `data` must be a data frame.", call. = FALSE)
  mode <- if (is.null(y)) "matrix" else "cross"
  x <- keep_existing(x, names(data), "cor_compute")
  y <- if (mode == "cross") keep_existing(y, names(data), "cor_compute") else NULL
  covariates <- if (length(covariates)) keep_existing(covariates, names(data), "cor_compute") else NULL
  partial <- length(covariates) > 0
  if (partial && !method %in% .core_methods)
    stop(sprintf("cor_compute: method \"%s\" has no partial version; use spearman, pearson or kendall.", method),
         call. = FALSE)
  if (partial && method == "kendall" && !requireNamespace("ppcor", quietly = TRUE))
    stop("cor_compute: partial Kendall correlations need the ppcor package.", call. = FALSE)
  if (!method %in% .core_methods && !requireNamespace("correlation", quietly = TRUE))
    stop(sprintf("cor_compute: method \"%s\" needs the correlation package.", method), call. = FALSE)

  # -- variables: type, numeric coding, usable? --------------------------------
  all_vars <- unique(c(x, y))
  types <- vapply(all_vars, function(v) cor_var_type(data[[v]]), character(1))
  num   <- lapply(stats::setNames(all_vars, all_vars), function(v) as_cor_numeric(data[[v]], v, method))
  usable <- vapply(num, function(v) sum(!is.na(v)) >= 5 && isTRUE(stats::sd(v, na.rm = TRUE) > 0), logical(1))
  if (any(!usable))
    warning("cor_compute: dropped (fewer than 5 values or constant): ",
            paste(all_vars[!usable], collapse = ", "), call. = FALSE)
  x <- x[usable[x]]
  if (mode == "cross") y <- y[usable[y]]
  if (mode == "matrix" && length(x) < 2) stop("cor_compute: fewer than 2 usable variables.", call. = FALSE)
  if (mode == "cross" && (!length(x) || !length(y))) stop("cor_compute: no usable x or y.", call. = FALSE)

  # -- pairs ----------------------------------------------------------------
  pairs <- if (mode == "matrix") {
    ij <- utils::combn(length(x), 2)
    tibble::tibble(X = x[ij[1, ]], Y = x[ij[2, ]])
  } else {
    tidyr::expand_grid(X = x, Y = y) %>% dplyr::filter(X != Y)
  }

  # -- computation ------------------------------------------------------------
  fast <- mode == "matrix" && !partial && engine %in% c("auto", "psych", "correlation", "stats")
  if (fast) {
    res <- .matrix_engine(as.data.frame(num[x], check.names = FALSE), method, engine, ci)
    tab <- dplyr::left_join(pairs, res, by = c("X", "Y"))
  } else {
    pair_engine <- if (engine %in% c("ppcor", "base")) engine else "ppcor"
    Z <- if (partial) data[, covariates, drop = FALSE] else NULL
    out <- purrr::map2(pairs$X, pairs$Y, function(a, b) {
      if (method %in% .core_methods) partial_cor_test(num[[a]], num[[b]], Z, method, pair_engine, ci)
      else .pair_extra(num[[a]], num[[b]], method, ci)
    })
    tab <- pairs %>%
      dplyr::mutate(r = purrr::map_dbl(out, "estimate"),
                    ci_low = purrr::map_dbl(out, "ci_low"), ci_high = purrr::map_dbl(out, "ci_high"),
                    p = purrr::map_dbl(out, "p.value"),
                    n = purrr::map_int(out, ~ as.integer(.x$n)),
                    k = purrr::map_int(out, ~ as.integer(.x$k)),
                    engine = purrr::map_chr(out, "engine"),
                    ci_engine = purrr::map_chr(out, "ci_engine"))
  }
  tab <- tab %>%
    dplyr::mutate(x_type = unname(types[X]), y_type = unname(types[Y]), .after = Y)

  # -- multiple-testing adjustment ---------------------------------------------
  adjusted <- !identical(p_adjust, "none")
  if (mode == "matrix") adjust_by <- "all"
  tab <- switch(if (adjusted) adjust_by else "none",
                y = dplyr::group_by(tab, Y), x = dplyr::group_by(tab, X), tab) %>%
    dplyr::mutate(p_adj = if (adjusted) stats::p.adjust(p, p_adjust) else p, .after = p) %>%
    dplyr::ungroup()
  family <- if (!adjusted) NULL else if (mode == "matrix") sprintf("the %d distinct pairs", nrow(tab)) else
    switch(adjust_by,
           y = sprintf("the %d predictors of each column", length(x)),
           x = sprintf("the %d columns of each predictor", length(y)),
           all = sprintf("all %d tests", nrow(tab)))

  meta <- list(mode = mode, meth = method, p_adjust = if (adjusted) p_adjust else "none",
               adjust_by = if (adjusted) adjust_by else "none", family = family,
               covariates = covariates, partial = partial,
               k = if (partial) suppressWarnings(max(tab$k, na.rm = TRUE)) else 0L,
               x_vars = x, y_vars = if (mode == "matrix") x else y,
               x_types = types[x], y_types = if (mode == "matrix") types[x] else types[y],
               labels = labels %||% attr(data, "labels"), bin_mark = bin_mark, ci = ci,
               engine = unique(tab$engine),
               ci_engine = setdiff(unique(tab$ci_engine), c("none", NA)))
  keep_cols <- unique(c(all_vars, covariates))
  structure(list(table = tab, meta = meta, data = data[, keep_cols, drop = FALSE]),
            class = c("cor_result", "list"))
}


# --------------------------------------------------------------------
# Internal engines
# --------------------------------------------------------------------

# Every pair of the columns of X (numeric data frame), raw p-values.
# Returns a tibble X, Y (both orders), r, ci_low, ci_high, p, n, k, engine, ci_engine.
.matrix_engine <- function(X, method, engine, ci) {
  vars <- names(X); k <- length(vars)
  has <- function(p) requireNamespace(p, quietly = TRUE)
  if (engine == "auto")
    engine <- if (method %in% c("pearson", "spearman") && has("psych")) "psych"
              else if (has("correlation")) "correlation" else "stats"
  if (engine == "psych" && method == "kendall") {
    message("cor_compute: psych tests Kendall's tau with the t formula of r; using stats::cor.test() instead.")
    engine <- "stats"
  }
  if (!method %in% .core_methods) engine <- "correlation"
  if (engine %in% c("psych", "correlation") && !has(engine))
    stop(sprintf("cor_compute: package %s not installed.", engine), call. = FALSE)

  R <- P <- N <- Lo <- Hi <- matrix(NA_real_, k, k, dimnames = list(vars, vars))
  if (engine == "psych") {
    ct <- psych::corr.test(X, use = "pairwise", method = method, adjust = "none", ci = FALSE)
    R <- ct$r; P <- ct$p
    N <- if (is.matrix(ct$n)) ct$n else matrix(ct$n, k, k, dimnames = list(vars, vars))
    eng <- "psych::corr.test"
  } else if (engine == "correlation") {
    cr <- as.data.frame(correlation::correlation(X, method = method, p_adjust = "none",
                                                 ci = ci %||% 0.95, redundant = TRUE))
    est <- intersect(c("r", "rho", "tau"), names(cr))[1]
    if (is.na(est)) stop("cor_compute: coefficient column not found in the output of correlation::correlation().",
                         call. = FALSE)
    i1 <- match(cr$Parameter1, vars); i2 <- match(cr$Parameter2, vars)
    ok <- !is.na(i1) & !is.na(i2) & i1 != i2
    idx <- cbind(i1[ok], i2[ok])
    R[idx] <- cr[[est]][ok]; P[idx] <- cr$p[ok]; N[idx] <- cr$n_Obs[ok]
    if (all(c("CI_low", "CI_high") %in% names(cr))) { Lo[idx] <- cr$CI_low[ok]; Hi[idx] <- cr$CI_high[ok] }
    eng <- "correlation::correlation"
  } else {
    for (i in seq_len(k)) for (j in seq_len(k)) {
      if (j <= i) next
      ok <- stats::complete.cases(X[[i]], X[[j]])
      if (sum(ok) < 5) next
      ct <- suppressWarnings(stats::cor.test(X[[i]][ok], X[[j]][ok], method = method, exact = FALSE))
      R[i, j] <- R[j, i] <- unname(ct$estimate)
      P[i, j] <- P[j, i] <- ct$p.value
      N[i, j] <- N[j, i] <- sum(ok)
    }
    eng <- "stats::cor.test"
  }
  ci_engine <- "none"
  if (!is.null(ci)) {
    if (eng == "correlation::correlation" && any(is.finite(Lo[upper.tri(Lo)]))) {
      ci_engine <- "correlation::correlation"
    } else if (method %in% .core_methods) {
      iv <- .cor_ci(as.vector(R), as.vector(N), method, ci)
      Lo[] <- iv$ci_low; Hi[] <- iv$ci_high
      ci_engine <- attr(iv, "engine")
    }
  }
  off <- row(R) != col(R)
  tibble::tibble(X = rownames(R)[row(R)[off]], Y = colnames(R)[col(R)[off]],
                 r = R[off], ci_low = Lo[off], ci_high = Hi[off], p = P[off],
                 n = as.integer(N[off]), k = 0L, engine = eng, ci_engine = ci_engine)
}

# One pair with a method of the correlation package (no partial version)
.pair_extra <- function(x, y, method, ci) {
  ok <- stats::complete.cases(x, y); n <- sum(ok)
  na <- list(estimate = NA_real_, p.value = NA_real_, n = n, k = 0L, ci_low = NA_real_, ci_high = NA_real_,
             engine = "none", ci_engine = "none")
  if (n < 5) return(na)
  out <- tryCatch(as.data.frame(correlation::cor_test(data.frame(x = x[ok], y = y[ok]), "x", "y",
                                                      method = method, ci = ci %||% 0.95)),
                  error = function(e) NULL)
  if (is.null(out)) return(na)
  est <- intersect(c("r", "rho", "tau"), names(out))[1]
  if (is.na(est)) return(na)
  has_ci <- !is.null(ci) && all(c("CI_low", "CI_high") %in% names(out))
  list(estimate = out[[est]][1], p.value = if ("p" %in% names(out)) out$p[1] else NA_real_, n = n, k = 0L,
       ci_low = if (has_ci) out$CI_low[1] else NA_real_, ci_high = if (has_ci) out$CI_high[1] else NA_real_,
       engine = "correlation::cor_test", ci_engine = if (has_ci) "correlation::correlation" else "none")
}


# --------------------------------------------------------------------
# Accessors and methods
# --------------------------------------------------------------------

# Internal: check the class
.check_result <- function(x, where) {
  if (!inherits(x, "cor_result")) stop(where, ": expected the output of cor_compute().", call. = FALSE)
  invisible(x)
}

#' Matrix of one statistic of a cor_result
#'
#' @param x output of [cor_compute()]
#' @param what column of `x$table`: `"r"` (default), `"p"`, `"p_adj"`, `"n"`,
#'   `"ci_low"` or `"ci_high"`
#' @param ... unused
#' @return a numeric matrix: symmetric with the diagonal filled (r = 1,
#'   p = 0, n = number of values) for a matrix analysis; `x` in rows and `y`
#'   in columns for a cross analysis
#' @export
as.matrix.cor_result <- function(x, what = c("r", "p", "p_adj", "n", "ci_low", "ci_high"), ...) {
  what <- match.arg(what)
  m <- x$meta; tab <- x$table
  rows <- m$x_vars; cols <- m$y_vars
  M <- matrix(NA_real_, length(rows), length(cols), dimnames = list(rows, cols))
  M[cbind(match(tab$X, rows), match(tab$Y, cols))] <- tab[[what]]
  if (m$mode == "matrix") {
    M[cbind(match(tab$Y, rows), match(tab$X, cols))] <- tab[[what]]
    diag(M) <- switch(what, r = 1, p = 0, p_adj = 0,
                      n = vapply(rows, function(v) sum(!is.na(x$data[[v]])), numeric(1)), NA_real_)
  }
  M
}

#' @export
as.data.frame.cor_result <- function(x, row.names = NULL, optional = FALSE, ...) as.data.frame(x$table)

#' @export
print.cor_result <- function(x, digits = 2, ...) {
  m <- x$meta; tab <- x$table
  n_bin <- sum(c(m$x_types, if (m$mode == "cross") m$y_types) == "binary")
  cat(sprintf("<cor_result> %s%s correlations: %s, %d pair%s%s\n",
              if (m$partial) "partial " else "", .meth_name(m$meth),
              if (m$mode == "matrix") sprintf("%d variables", length(m$x_vars))
              else sprintf("%d x %d variables", length(m$x_vars), length(m$y_vars)),
              nrow(tab), if (nrow(tab) > 1) "s" else "",
              if (n_bin) sprintf(" (%d binary variable%s)", n_bin, if (n_bin > 1) "s" else "") else ""))
  if (m$partial) cat("  adjusted for:", paste(m$covariates, collapse = ", "), "\n")
  cat(sprintf("  engine: %s | N = %s\n", .cor_engine_txt(tab$engine), .n_range_txt(tab$n)))
  cat("  ", .cor_p_line(m$p_adjust, m$family), "\n\n", sep = "")
  R <- as.matrix(x, "r"); P <- as.matrix(x, "p_adj")
  lab <- matrix(ifelse(is.na(R), "", paste0(sprintf(paste0("%.", digits, "f"), R), p_to_stars(P, ns = ""))),
                nrow(R), dimnames = dimnames(R))
  if (m$mode == "matrix") { lab[upper.tri(lab)] <- ""; diag(lab) <- "1" }
  print(noquote(lab))
  invisible(x)
}

# Internal: subtitle shared by the plots and the tables
.cor_subtitle <- function(meta, tab, labels = meta$labels, newline = "\n") {
  sprintf("%s correlation, %s | %sN = %s%sComputed with %s; p adjustment: %s",
          .meth_name(meta$meth),
          if (meta$partial) paste0("partial, adjusted for ", paste(get_label(meta$covariates, labels), collapse = ", "))
          else "unadjusted",
          if (meta$mode == "matrix") "pairwise complete " else "",
          .n_range_txt(tab$n), newline, .cor_engine_txt(tab$engine),
          if (identical(meta$p_adjust, "none")) "none (raw p-values)" else "stats::p.adjust() {stats}")
}
