# ====================================================================
# cor_pair.R -- one pair: simple or partial correlation test, with CI
# ====================================================================

#' Confidence interval of a correlation coefficient (Fisher z)
#'
#' Uses `correlation::cor_to_ci()` when the package is installed, otherwise the
#' same Fisher z transformation in base R: z = atanh(r), z +/- q SE, back with
#' tanh(), with SE = 1 / sqrt(n - 3) (Pearson), sqrt(1.06 / (n - 3)) (Spearman)
#' and sqrt(0.437 / (n - 4)) (Kendall) (Fieller, Hartley & Pearson 1957).
#' For a partial correlation pass `n - k` (k = number of covariate terms).
#' @return data.frame(ci_low, ci_high), attribute "engine"
#' @noRd
.cor_ci <- function(r, n, meth = "pearson", level = 0.95) {
  r <- as.numeric(r); n <- rep_len(as.numeric(n), length(r))
  out <- data.frame(ci_low = rep(NA_real_, length(r)), ci_high = rep(NA_real_, length(r)))
  attr(out, "engine") <- "none"
  if (is.null(level) || is.na(level) || !length(r)) return(out)
  ok <- is.finite(r) & abs(r) < 1 & is.finite(n) & n > 4
  if (!any(ok)) return(out)
  if (meth %in% .core_methods && requireNamespace("correlation", quietly = TRUE)) {
    res <- tryCatch(correlation::cor_to_ci(r[ok], n = n[ok], ci = level, method = meth), error = function(e) NULL)
    if (!is.null(res) && length(res$CI_low) == sum(ok)) {
      out$ci_low[ok] <- res$CI_low; out$ci_high[ok] <- res$CI_high
      attr(out, "engine") <- "correlation::cor_to_ci"
      return(out)
    }
  }
  se <- switch(meth, spearman = sqrt(1.06 / (n[ok] - 3)), kendall = sqrt(0.437 / (n[ok] - 4)), 1 / sqrt(n[ok] - 3))
  q  <- stats::qnorm(1 - (1 - level) / 2)
  z  <- atanh(r[ok])
  out$ci_low[ok] <- tanh(z - q * se); out$ci_high[ok] <- tanh(z + q * se)
  attr(out, "engine") <- "base R"
  out
}

#' Rows used and numeric design of the covariates
#'
#' Complete cases of (x, y, Z); covariates converted with model.matrix()
#' (factors -> dummies); covariates with a single value and constant columns dropped.
#' @return list(ok = logical vector, Zm = matrix or NULL, k = number of columns of Zm)
#' @noRd
.partial_design <- function(x, y, Z = NULL) {
  Zd <- if (!is.null(Z)) as.data.frame(Z) else NULL
  if (!is.null(Zd) && ncol(Zd) == 0) Zd <- NULL
  ok <- if (is.null(Zd)) stats::complete.cases(x, y) else stats::complete.cases(x, y, Zd)
  Zm <- NULL
  if (!is.null(Zd) && sum(ok) > 0) {
    Zs <- droplevels(Zd[ok, , drop = FALSE])
    Zs <- Zs[, vapply(Zs, function(v) length(unique(v)) > 1, logical(1)), drop = FALSE]
    if (ncol(Zs) > 0) {
      Zm <- stats::model.matrix(~ ., data = Zs)[, -1, drop = FALSE]
      Zm <- Zm[, apply(Zm, 2, stats::sd) > 0, drop = FALSE]
      if (ncol(Zm) == 0) Zm <- NULL
    }
  }
  list(ok = ok, Zm = Zm, k = if (is.null(Zm)) 0L else ncol(Zm))
}

#' Residuals of x and y on the covariates (ranks everywhere for Spearman, as in ppcor)
#' @return list(ex, ey)
#' @noRd
.partial_residuals <- function(x, y, Zm, meth) {
  rk <- function(v) if (meth == "spearman") rank(v) else v
  Zf <- if (meth == "spearman") apply(Zm, 2, rank) else Zm
  D  <- cbind(1, Zf)
  list(ex = stats::lm.fit(D, rk(x))$residuals, ey = stats::lm.fit(D, rk(y))$residuals)
}

#' Correlation test between x and y, partial if covariates are given
#'
#' * no covariate: [stats::cor.test()] (`exact = FALSE`);
#' * covariates: `ppcor::pcor.test()` (default engine); fallback or
#'   `engine = "base"`: correlation of the residuals of x and y on Z
#'   (for Spearman x, y AND Z are ranked, as in ppcor), t-test with
#'   n - 2 - k df. For Pearson this is exactly the slope test of
#'   `lm(y ~ x + Z)` (Frisch-Waugh-Lovell theorem). Kendall partial
#'   correlations need ppcor.
#'
#' @param x,y numeric vectors (binary: 0/1, see [as_cor_numeric()])
#' @param Z data frame of covariates (factors allowed: dummy coded) or `NULL`
#' @param meth `"spearman"` (ranks), `"pearson"` (values) or `"kendall"`
#' @param engine `"ppcor"` (default) or `"base"` (residual method)
#' @param ci confidence level of the interval (`NULL` = none); Fisher z on n - k
#' @return list(estimate, p.value, n, k, ci_low, ci_high, engine, ci_engine)
#' @export
#' @examples
#' d <- simulate_cor_data(100)
#' partial_cor_test(d$IL6, d$fatigue, meth = "spearman")
#' partial_cor_test(d$NfL, d$cognition, Z = d[, c("age", "sex")], meth = "pearson")
partial_cor_test <- function(x, y, Z = NULL, meth = c("spearman", "pearson", "kendall"),
                             engine = c("ppcor", "base"), ci = 0.95) {
  meth <- match.arg(meth); engine <- match.arg(engine)
  des <- .partial_design(x, y, Z)
  x <- x[des$ok]; y <- y[des$ok]; n <- length(x)
  out <- function(r, p, k, eng) {
    iv <- .cor_ci(r, n - k, meth, ci)
    list(estimate = unname(r), p.value = unname(p), n = n, k = as.integer(k),
         ci_low = iv$ci_low, ci_high = iv$ci_high, engine = eng, ci_engine = attr(iv, "engine"))
  }
  if (n < 5 || stats::sd(x) == 0 || stats::sd(y) == 0) return(out(NA_real_, NA_real_, 0, "none"))

  simple <- function() {
    ct <- suppressWarnings(stats::cor.test(x, y, method = meth, exact = FALSE))
    out(ct$estimate, ct$p.value, 0, "stats::cor.test")
  }
  Zm <- des$Zm
  if (is.null(Zm)) return(simple())
  k <- ncol(Zm)
  if (n - 2 - k < 1) return(out(NA_real_, NA_real_, k, "none"))

  if (engine == "ppcor" && requireNamespace("ppcor", quietly = TRUE)) {
    res <- tryCatch(suppressWarnings(ppcor::pcor.test(x, y, Zm, method = meth)), error = function(e) NULL)
    if (!is.null(res) && is.finite(res$estimate)) return(out(res$estimate, res$p.value, k, "ppcor::pcor.test"))
  }
  if (meth == "kendall")
    stop("partial_cor_test: a partial Kendall correlation needs the ppcor package (install.packages(\"ppcor\")).",
         call. = FALSE)
  # Fallback (or engine = "base"): correlation of the residuals
  e  <- .partial_residuals(x, y, Zm, meth)
  r  <- stats::cor(e$ex, e$ey)
  df <- n - 2 - k
  t  <- r * sqrt(df / max(1e-12, 1 - r^2))
  out(r, 2 * stats::pt(-abs(t), df), k, "residuals (base R)")
}

#' @rdname partial_cor_test
#' @description `partial_rank_cor()` is the former name, kept for backward compatibility.
#' @export
partial_rank_cor <- function(x, y, Z = NULL, meth = c("spearman", "pearson"))
  partial_cor_test(x, y, Z, match.arg(meth))
