# ====================================================================
# binary_effects.R -- effect sizes to report for a binary (0/1) variable
# ====================================================================

#' Effect sizes of binary variables on continuous scores
#'
#' For a binary variable (e.g. a biomarker "high" vs "low"), the correlation is
#' bounded by the split between groups and is not comparable with the
#' correlation of a continuous variable (see [cor_binary_note()]). This
#' function returns what to report instead, with existing tests:
#' * difference in means M_1 - M_0 and its CI: [stats::t.test()] with equal
#'   variances (the test of the point-biserial r);
#' * rank-biserial correlation r_rb = 2U / (n_0 n_1) - 1 = P(Y_1 > Y_0) -
#'   P(Y_1 < Y_0) = 2 AUC - 1 (Kerby 2014), with U the Mann-Whitney statistic of
#'   [stats::wilcox.test()]; CI from `effectsize::rank_biserial()` if installed;
#' * for comparison, the observed r_s and r_pb and the bound sqrt(3pq) of |r_s|.
#'
#' Group 1 = upper level of the binary variable (as in [as_cor_numeric()]). p-values are raw.
#'
#' @param data data frame
#' @param bin_vars binary variables (non-binary ones are dropped with a warning)
#' @param y_vars continuous outcomes
#' @param labels named vector variable -> display label (kept for the tables)
#' @param ci confidence level
#' @return tibble, one row per (binary variable, score): X, Y, level0, level1,
#'   n0, n1, prop1, mean0, sd0, mean1, sd1, median0, median1, diff, diff_low,
#'   diff_high, p_t, r_rb, r_rb_low, r_rb_high, p_mw, r_s, r_s_bound, r_pb.
#' @references Kerby DS (2014). Comprehensive Psychology 3:11.IT.3.1. doi:10.2466/11.IT.3.1
#' @export
#' @examples
#' d <- simulate_cor_data(150)
#' cor_binary_effects(d, c("BIN_CRP", "BIN_IL6"), c("depression", "fatigue"))
cor_binary_effects <- function(data, bin_vars, y_vars, labels = NULL, ci = 0.95) {
  bin_vars <- keep_existing(bin_vars, names(data), "cor_binary_effects")
  y_vars   <- keep_existing(y_vars, names(data), "cor_binary_effects")
  types <- vapply(bin_vars, function(v) cor_var_type(data[[v]]), character(1))
  if (any(types != "binary"))
    warning("cor_binary_effects: not binary, ignored: ", paste(bin_vars[types != "binary"], collapse = ", "),
            call. = FALSE)
  bin_vars <- bin_vars[types == "binary"]
  if (!length(bin_vars) || !length(y_vars)) stop("cor_binary_effects: no usable variable.", call. = FALSE)

  rows <- list()
  for (b in bin_vars) {
    lev <- .binary_levels(data[[b]])
    gv  <- as_cor_numeric(data[[b]], b, "spearman")
    for (yv in y_vars) {
      yy <- as.numeric(data[[yv]])
      ok <- stats::complete.cases(gv, yy)
      y1 <- yy[ok & gv == 1]; y0 <- yy[ok & gv == 0]
      n1 <- length(y1); n0 <- length(y0)
      row <- tibble::tibble(X = b, Y = yv, level0 = lev[1], level1 = lev[2], n0 = n0, n1 = n1,
                            prop1 = n1 / (n0 + n1),
                            mean0 = mean(y0), sd0 = stats::sd(y0), mean1 = mean(y1), sd1 = stats::sd(y1),
                            median0 = stats::median(y0), median1 = stats::median(y1),
                            diff = NA_real_, diff_low = NA_real_, diff_high = NA_real_, p_t = NA_real_,
                            r_rb = NA_real_, r_rb_low = NA_real_, r_rb_high = NA_real_, p_mw = NA_real_,
                            r_s = NA_real_, r_s_bound = NA_real_, r_pb = NA_real_)
      if (n1 >= 2 && n0 >= 2 && stats::sd(yy[ok]) > 0) {
        tt <- stats::t.test(y1, y0, var.equal = TRUE, conf.level = ci)
        wt <- suppressWarnings(stats::wilcox.test(y1, y0, exact = FALSE))
        U  <- unname(wt$statistic)
        rb <- 2 * U / (n1 * n0) - 1
        rb_ci <- .rb_ci(y1, y0, ci, rb)
        p1 <- n1 / (n0 + n1)
        row$diff <- mean(y1) - mean(y0)
        row$diff_low <- tt$conf.int[1]; row$diff_high <- tt$conf.int[2]; row$p_t <- tt$p.value
        row$r_rb <- rb; row$r_rb_low <- rb_ci[1]; row$r_rb_high <- rb_ci[2]; row$p_mw <- wt$p.value
        row$r_s <- stats::cor(gv[ok], yy[ok], method = "spearman")
        row$r_pb <- stats::cor(gv[ok], yy[ok])
        row$r_s_bound <- sqrt(3 * p1 * (1 - p1))
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  out <- dplyr::bind_rows(rows)
  attr(out, "cor_meta") <- list(labels = labels, ci = ci,
                                rb_ci_engine = if (requireNamespace("effectsize", quietly = TRUE))
                                  "effectsize::rank_biserial" else "none")
  out
}

# Internal: CI of r_rb from effectsize::rank_biserial(); its sign is aligned on r_rb (group 1 - group 0)
.rb_ci <- function(y1, y0, ci, r_rb) {
  if (!requireNamespace("effectsize", quietly = TRUE)) return(c(NA_real_, NA_real_))
  rb <- tryCatch(as.data.frame(effectsize::rank_biserial(y1, y0, ci = ci)), error = function(e) NULL)
  if (is.null(rb) || !all(c("CI_low", "CI_high") %in% names(rb))) return(c(NA_real_, NA_real_))
  lo <- rb$CI_low[1]; hi <- rb$CI_high[1]
  est <- if ("r_rank_biserial" %in% names(rb)) rb$r_rank_biserial[1] else NA_real_
  if (is.finite(est) && is.finite(r_rb) && abs(r_rb) > 1e-12 && sign(est) != sign(r_rb)) {
    tmp <- -hi; hi <- -lo; lo <- tmp
  }
  c(lo, hi)
}
