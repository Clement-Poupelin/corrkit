# ====================================================================
# plot_pair.R -- one pair: what the coefficient summarises
# ====================================================================

#' Plot of one correlation, with its coefficient, CI and test in the subtitle
#'
#' What is drawn is what the coefficient measures:
#' * two continuous variables, unadjusted: scatter plot of the raw values,
#'   least-squares line (Pearson) or LOESS curve (Spearman, Kendall: monotonic trend);
#' * one binary variable, unadjusted: the continuous variable by group (box,
#'   points, diamond = mean), with r_rb and the difference in means
#'   ([cor_binary_effects()]);
#' * with covariates: added-variable plot, i.e. residuals of y against
#'   residuals of x on the covariates (ranks for Spearman), whose correlation
#'   IS the partial correlation.
#'
#' @inheritParams partial_cor_test
#' @param data data frame
#' @param x,y names of the two variables
#' @param covariates covariates of the partial correlation (`NULL` = unadjusted)
#' @param labels named vector variable -> display label
#' @param bin_mark marker appended to the label of a binary variable
#' @param base_size base font size
#' @param title title (default "<y> vs <x>")
#' @param palette colours (named low / high); `high` is used for the line and the means
#' @return a ggplot
#' @export
#' @examples
#' d <- simulate_cor_data(150)
#' plot_cor_pair(d, "IL6", "fatigue", labels = attr(d, "labels"))
#' plot_cor_pair(d, "BIN_CRP", "depression", labels = attr(d, "labels"))
#' plot_cor_pair(d, "NfL", "cognition", covariates = "age", meth = "pearson")
plot_cor_pair <- function(data, x, y, covariates = NULL, labels = NULL,
                          meth = c("spearman", "pearson", "kendall"), engine = c("ppcor", "base"),
                          ci = 0.95, bin_mark = " [0/1]", base_size = 11, title = NULL, palette = COR_PAL) {
  meth <- match.arg(meth); engine <- match.arg(engine)
  if (!all(c(x, y) %in% names(data))) stop("plot_cor_pair: x or y not found in the data.", call. = FALSE)
  covariates <- if (length(covariates)) keep_existing(covariates, names(data), "plot_cor_pair") else NULL
  Z  <- if (length(covariates)) data[, covariates, drop = FALSE] else NULL
  xt <- cor_var_type(data[[x]]); yt <- cor_var_type(data[[y]])
  xn <- as_cor_numeric(data[[x]], x, meth); yn <- as_cor_numeric(data[[y]], y, meth)
  res <- partial_cor_test(xn, yn, Z, meth, engine, ci)
  des <- .partial_design(xn, yn, Z)
  partial <- des$k > 0
  xl <- mark_binary(get_label(x, labels), xt, bin_mark)
  yl <- mark_binary(get_label(y, labels), yt, bin_mark)
  sym <- .cor_sym_txt(meth)
  ci_txt <- if (!is.null(ci) && is.finite(res$ci_low))
    sprintf(" [%.2f, %.2f]", res$ci_low, res$ci_high) else ""
  stat_txt <- sprintf("%s%s = %.2f%s, %s, n = %d%s", if (partial) "partial " else "", sym, res$estimate, ci_txt,
                      .fmt_p_eq(res$p.value), res$n,
                      if (partial) paste0(" | adjusted for ", paste(get_label(covariates, labels), collapse = ", ")) else "")
  eng_txt <- sprintf("Computed with %s.", .cor_engine_txt(res$engine))
  ci_lvl  <- if (!is.null(ci)) sprintf("%g%% CI", 100 * ci) else "CI"
  thm <- theme_pub(base_size) + .cor_text_theme()

  # (a) partial: added-variable plot -----------------------------------------
  if (partial) {
    ok <- des$ok
    e  <- .partial_residuals(xn[ok], yn[ok], des$Zm, meth)
    d  <- data.frame(ex = e$ex, ey = e$ey)
    rank_txt <- if (meth == "spearman") " (ranks)" else ""
    return(ggplot2::ggplot(d, ggplot2::aes(x = ex, y = ey)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey80") + ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
      ggplot2::geom_point(alpha = 0.6, size = 1.8, colour = "grey25") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = palette[["high"]],
                           fill = palette[["high"]], alpha = 0.15) +
      ggplot2::labs(title = title %||% sprintf("%s vs %s, adjusted", yl, xl),
                    subtitle = paste0(stat_txt, "\n", eng_txt),
                    x = sprintf("%s%s, residual on the covariates", xl, rank_txt),
                    y = sprintf("%s%s, residual on the covariates", yl, rank_txt),
                    caption = .wrap_caption(list(
                      paste("How to read: each point is an individual once the covariates have been removed from both variables",
                            "(residuals of a linear regression on the covariates). The correlation of these residuals is the partial",
                            "correlation; line = least-squares fit with its 95% confidence band."),
                      .cor_formula_lines(meth, TRUE, des$k)), 110)) +
      thm)
  }

  # (b) one binary variable: values by group ---------------------------------
  if (xt == "binary" || yt == "binary") {
    b <- if (xt == "binary") x else y
    v <- if (xt == "binary") y else x
    bl <- if (xt == "binary") xl else yl
    vl <- get_label(v, labels)
    eff <- cor_binary_effects(data, b, v, ci = ci %||% 0.95)
    lev <- c(eff$level0, eff$level1)
    gn <- as_cor_numeric(data[[b]], b, meth); vv <- as.numeric(data[[v]])
    ok <- stats::complete.cases(gn, vv)
    d <- data.frame(g = factor(ifelse(gn[ok] == 1, lev[2], lev[1]), levels = lev), v = vv[ok])
    rb_ci <- if (is.finite(eff$r_rb_low)) sprintf(" [%.2f, %.2f]", eff$r_rb_low, eff$r_rb_high) else ""
    eff_txt <- sprintf("r_rb = %.2f%s | mean difference %s \u2212 %s = %.2f [%.2f, %.2f] | max |r_s| for this split = %.2f",
                       eff$r_rb, rb_ci, lev[2], lev[1], eff$diff, eff$diff_low, eff$diff_high, eff$r_s_bound)
    return(ggplot2::ggplot(d, ggplot2::aes(x = g, y = v)) +
      ggplot2::geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey93", colour = "grey40") +
      ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.5, size = 1.6, colour = "grey25") +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = palette[["high"]], colour = "white") +
      ggplot2::scale_x_discrete(labels = sprintf("%s\n(n = %d)", lev, c(eff$n0, eff$n1))) +
      ggplot2::labs(title = title %||% sprintf("%s by %s", vl, bl),
                    subtitle = paste0(stat_txt, "\n", eff_txt, "\n", eng_txt),
                    x = bl, y = vl,
                    caption = .wrap_caption(list(
                      c("How to read: points = individuals; box = median and quartiles; diamond = mean of the group.",
                        sprintf("r_rb = P(%s > %s) \u2212 P(%s < %s) for a random pair of individuals (= 2 AUC \u2212 1); mean difference with its %s (Student's t-test).",
                                lev[2], lev[1], lev[2], lev[1], ci_lvl)),
                      cor_binary_note(meth, "en", short = TRUE, mark = bin_mark)), 110)) +
      thm)
  }

  # (c) two continuous variables ---------------------------------------------
  ok <- stats::complete.cases(xn, yn)
  d  <- data.frame(ex = xn[ok], ey = yn[ok])
  smooth <- if (meth == "pearson") "lm" else "loess"
  ggplot2::ggplot(d, ggplot2::aes(x = ex, y = ey)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8, colour = "grey25") +
    ggplot2::geom_smooth(method = smooth, formula = y ~ x, se = TRUE, colour = palette[["high"]],
                         fill = palette[["high"]], alpha = 0.15) +
    ggplot2::labs(title = title %||% sprintf("%s vs %s", yl, xl),
                  subtitle = paste0(stat_txt, "\n", eng_txt), x = xl, y = yl,
                  caption = .wrap_caption(list(
                    c("How to read: each point is an individual.",
                      if (meth == "pearson") "Line = least-squares fit with its 95% confidence band: r measures how close the points are to a straight line."
                      else sprintf("Curve = LOESS smoother with its 95%% band: %s measures a monotonic (not necessarily linear) trend, computed on the ranks.", sym)),
                    .cor_formula_lines(meth)), 110)) +
    thm
}
