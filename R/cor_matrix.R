# ====================================================================
# cor_matrix.R -- matrix of all pairs
# ====================================================================

#' Correlation and p-value matrices of all pairs (pairwise complete observations)
#'
#' The engine always returns **raw** p-values. By default they are kept as
#' they are (`p_adjust = "none"`); with another `p_adjust`, corrkit adjusts
#' them with [stats::p.adjust()] over the k(k - 1)/2 distinct pairs, as
#' `psych::corr.test()` does above its diagonal.
#'
#' Engines (`engine = "auto"` picks the first available in this order):
#' * Pearson / Spearman: `psych::corr.test(adjust = "none")`, then
#'   `correlation::correlation(p_adjust = "none")`, then [stats::cor.test()] pair by pair;
#' * Kendall: `correlation::correlation()`, then [stats::cor.test()] (psych
#'   tests tau with the t formula of r, so it is not used for Kendall);
#' * other methods of the correlation package (`"biweight"`, `"percentage"`,
#'   `"blomqvist"`, `"hoeffding"`, `"gaussian"`, `"shepherd"`, `"distance"`):
#'   `correlation::correlation()` only.
#'
#' Binary variables are coded 0/1 (see [as_cor_numeric()]) and flagged in `types`.
#'
#' @param data data frame
#' @param vars variables to correlate
#' @param method correlation method (see above)
#' @param p_adjust method of [stats::p.adjust()]: `"none"` (default, raw p-values, as in most
#'   correlation analyses), or e.g. `"BH"`, `"holm"`, `"bonferroni"`
#' @param engine `"auto"`, `"psych"`, `"correlation"` or `"stats"`
#' @param ci confidence level of the intervals (`NULL` = none)
#' @return object of class `cor_matrix`: list(r, p, p_adj, n, ci_low, ci_high,
#'   long, types, method, p_adjust, engine, ci, ci_engine); `long` = one row per
#'   pair (i < j) with Var1, Var2, type1, type2, r, ci_low, ci_high, p, p_adj, n.
#' @export
#' @examples
#' d  <- simulate_cor_data(120)
#' cm <- cor_test_matrix(d, c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP"), method = "spearman")
#' cm
#' head(cm$long)
cor_test_matrix <- function(data, vars, method = "spearman", p_adjust = "none",
                            engine = c("auto", "psych", "correlation", "stats"), ci = 0.95) {
  engine <- match.arg(engine)
  method <- match.arg(method, c(.core_methods, .extra_methods))
  vars <- keep_existing(vars, names(data), "cor_test_matrix")
  types <- vapply(vars, function(v) cor_var_type(data[[v]]), character(1))
  X <- as.data.frame(lapply(stats::setNames(vars, vars), function(v) as_cor_numeric(data[[v]], v, method)),
                     check.names = FALSE)
  keep <- vapply(X, function(v) sum(!is.na(v)) >= 5 && stats::sd(v, na.rm = TRUE) > 0, logical(1))
  if (any(!keep))
    warning("cor_test_matrix: dropped (fewer than 5 values or constant): ", paste(names(X)[!keep], collapse = ", "),
            call. = FALSE)
  X <- X[, keep, drop = FALSE]; types <- types[keep]
  vars <- names(X); k <- length(vars)
  if (k < 2) stop("cor_test_matrix: fewer than 2 usable variables.", call. = FALSE)

  has <- function(p) requireNamespace(p, quietly = TRUE)
  if (engine == "auto") {
    engine <- if (method %in% c("pearson", "spearman") && has("psych")) "psych"
              else if (has("correlation")) "correlation" else "stats"
  }
  if (engine == "psych" && method == "kendall") {
    message("cor_test_matrix: psych tests Kendall's tau with the t formula of r; using stats::cor.test() instead.")
    engine <- "stats"
  }
  if (!method %in% .core_methods && engine != "correlation")
    stop(sprintf("cor_test_matrix: method \"%s\" needs engine = \"correlation\" (package correlation).", method),
         call. = FALSE)
  if (engine %in% c("psych", "correlation") && !has(engine))
    stop(sprintf("cor_test_matrix: package %s not installed.", engine), call. = FALSE)

  empty <- function() matrix(NA_real_, k, k, dimnames = list(vars, vars))
  R <- P <- N <- Lo <- Hi <- empty()

  if (engine == "psych") {
    ct <- psych::corr.test(X, use = "pairwise", method = method, adjust = "none", ci = FALSE)
    R <- ct$r; P <- ct$p
    N <- if (is.matrix(ct$n)) ct$n else matrix(ct$n, k, k, dimnames = list(vars, vars))
    eng_name <- "psych::corr.test"
  } else if (engine == "correlation") {
    cr <- as.data.frame(correlation::correlation(X, method = method, p_adjust = "none",
                                                 ci = ci %||% 0.95, redundant = TRUE))
    est <- intersect(c("r", "rho", "tau"), names(cr))[1]
    if (is.na(est)) stop("cor_test_matrix: coefficient column not found in the output of correlation::correlation().",
                         call. = FALSE)
    i1 <- match(cr$Parameter1, vars); i2 <- match(cr$Parameter2, vars)
    ok <- !is.na(i1) & !is.na(i2) & i1 != i2
    idx <- cbind(i1[ok], i2[ok])
    R[idx] <- cr[[est]][ok]; P[idx] <- cr$p[ok]; N[idx] <- cr$n_Obs[ok]
    if (all(c("CI_low", "CI_high") %in% names(cr))) { Lo[idx] <- cr$CI_low[ok]; Hi[idx] <- cr$CI_high[ok] }
    eng_name <- "correlation::correlation"
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
    eng_name <- "stats::cor.test"
  }
  diag(R) <- 1; diag(P) <- 0
  diag(N) <- vapply(X, function(v) sum(!is.na(v)), numeric(1))

  # Confidence intervals: from correlation::correlation() when it computed them, else Fisher z
  ci_engine <- "none"
  if (!is.null(ci)) {
    if (eng_name == "correlation::correlation" && any(is.finite(Lo[upper.tri(Lo)]))) {
      ci_engine <- "correlation::correlation"
    } else if (method %in% .core_methods) {
      up <- upper.tri(R)
      iv <- .cor_ci(R[up], N[up], method, ci)
      Lo[up] <- iv$ci_low; Hi[up] <- iv$ci_high
      Lo[lower.tri(Lo)] <- t(Lo)[lower.tri(Lo)]; Hi[lower.tri(Hi)] <- t(Hi)[lower.tri(Hi)]
      ci_engine <- attr(iv, "engine")
    }
  }

  up <- upper.tri(P)
  long <- tibble::tibble(Var1 = rownames(P)[row(P)[up]], Var2 = colnames(P)[col(P)[up]],
                         type1 = unname(types[Var1]), type2 = unname(types[Var2]),
                         r = R[up], ci_low = Lo[up], ci_high = Hi[up], p = P[up], n = N[up]) %>%
    dplyr::mutate(p_adj = stats::p.adjust(p, p_adjust))
  Padj <- P
  Padj[up] <- long$p_adj
  Padj[lower.tri(Padj)] <- t(Padj)[lower.tri(Padj)]
  structure(list(r = R, p = P, p_adj = Padj, n = N, ci_low = Lo, ci_high = Hi, long = long, types = types,
                 method = method, p_adjust = p_adjust, engine = eng_name, ci = ci, ci_engine = ci_engine),
            class = c("cor_matrix", "list"))
}

#' @export
print.cor_matrix <- function(x, digits = 2, ...) {
  k <- nrow(x$r)
  cat(sprintf("<cor_matrix> %s correlations, %d variables (%d binary), %d pairs\n",
              .meth_name(x$method), k, sum(x$types == "binary"), k * (k - 1) / 2))
  cat(sprintf("  engine: %s | %s | N = %s\n", .cor_engine_txt(x$engine),
              if (identical(x$p_adjust, "none")) "raw p-values (no adjustment)"
              else sprintf("p adjustment: %s over the %d pairs", .p_adjust_name(x$p_adjust), k * (k - 1) / 2),
              .n_range_txt(x$n[upper.tri(x$n)])))
  lab <- matrix(paste0(sprintf(paste0("%.", digits, "f"), x$r), p_to_stars(x$p_adj, ns = "")), k, k,
                dimnames = dimnames(x$r))
  lab[upper.tri(lab)] <- ""
  diag(lab) <- "1"
  print(noquote(lab))
  invisible(x)
}
