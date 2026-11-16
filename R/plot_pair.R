# ====================================================================
# plot_pair.R -- one pair: what the coefficient summarises
# ====================================================================

#' Plot of one correlation, with its coefficient, CI and test in the subtitle
#'
#' What is drawn is what the coefficient measures:
#' * two continuous variables, unadjusted: scatter plot of the raw values with
#'   a fitted curve; by default a least-squares line for Pearson (r measures a
#'   straight line) and a LOESS curve for Spearman and Kendall (monotonic trend);
#' * one binary variable, unadjusted: the continuous variable by group (box,
#'   points, diamond = mean), with r_rb and the difference in means
#'   ([cor_binary_effects()]);
#' * with covariates: added-variable plot, i.e. residuals of y against
#'   residuals of x on the covariates (ranks for Spearman), whose correlation
#'   IS the partial correlation; by default a least-squares line.
#'
#' Default colours come from the Okabe-Ito palette, readable with the common
#' colour-vision deficiencies.
#'
#' @inheritParams partial_cor_test
#' @param data data frame
#' @param x,y names of the two variables
#' @param covariates covariates of the partial correlation (`NULL` = unadjusted)
#' @param labels named vector variable -> display label
#' @param smooth fitted curve, passed to `ggplot2::geom_smooth(method = )`:
#'   `"lm"`, `"loess"`, `"gam"`, `"glm"`, a function, or `"none"` for no curve.
#'   `NULL` (default) = `"lm"` for Pearson and partial correlations, `"loess"`
#'   for Spearman and Kendall.
#' @param smooth_args list of other arguments of `ggplot2::geom_smooth()`,
#'   e.g. `list(span = 0.5)` for LOESS or `list(formula = y ~ poly(x, 2))`
#' @param se draw the confidence band of the curve
#' @param point_colour colour of the points (scatter plots)
#' @param line_colour colour of the curve and its band
#' @param group_colours colours of the two groups (binary variable), lower level first
#' @param bin_mark marker appended to the label of a binary variable
#' @param base_size base font size
#' @param title title (default "<y> vs <x>")
#' @return a ggplot
#' @export
#' @examples
#' d <- simulate_cor_data(150)
#' plot_cor_pair(d, "IL6", "fatigue", labels = attr(d, "labels"))
#' plot_cor_pair(d, "IL6", "fatigue", smooth = "lm")
#' plot_cor_pair(d, "IL6", "fatigue", smooth = "loess", smooth_args = list(span = 0.5))
#' plot_cor_pair(d, "BIN_CRP", "depression", labels = attr(d, "labels"))
#' plot_cor_pair(d, "NfL", "cognition", covariates = "age", meth = "pearson")
plot_cor_pair <- function(data, x, y, covariates = NULL, labels = NULL,
                          meth = c("spearman", "pearson", "kendall"), engine = c("ppcor", "base"),
                          ci = 0.95, smooth = NULL, smooth_args = list(), se = TRUE,
                          point_colour = "#0072B2", line_colour = "#D55E00",
                          group_colours = c("#56B4E9", "#E69F00"),
                          bin_mark = " [0/1]", base_size = 11, title = NULL) {
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
  cov_txt <- paste(get_label(covariates, labels), collapse = ", ")
  sym <- .cor_sym_txt(meth)
  ci_txt <- if (!is.null(ci) && is.finite(res$ci_low))
    sprintf(" [%.2f, %.2f]", res$ci_low, res$ci_high) else ""
  stat_txt <- sprintf("%s%s = %.2f%s, %s, n = %d%s", if (partial) "partial " else "", sym, res$estimate, ci_txt,
                      .fmt_p_eq(res$p.value), res$n, if (partial) paste0(" | adjusted for ", cov_txt) else "")
  eng_txt <- sprintf("Computed with %s.", .cor_engine_txt(res$engine))
  ci_lvl  <- if (!is.null(ci)) sprintf("%g%% CI", 100 * ci) else "CI"
  thm <- theme_pub(base_size) + .cor_text_theme()
  smooth <- smooth %||% (if (partial || meth == "pearson") "lm" else "loess")
  fit <- .smooth_layer(smooth, smooth_args, line_colour, se)

  # (a) partial: added-variable plot -----------------------------------------
  if (partial) {
    ok <- des$ok
    e  <- .partial_residuals(xn[ok], yn[ok], des$Zm, meth)
    d  <- data.frame(ex = e$ex, ey = e$ey)
    what <- if (meth == "spearman") "rank residual" else "residual"
    return(ggplot2::ggplot(d, ggplot2::aes(x = ex, y = ey)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey80") + ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
      ggplot2::geom_point(alpha = 0.6, size = 1.8, colour = point_colour) +
      fit$layer +
      .labs_safe(title = title %||% sprintf("%s vs %s, adjusted", yl, xl),
                 subtitle = paste0(stat_txt, "\n", eng_txt),
                 x = sprintf("%s (%s)", xl, what), y = sprintf("%s (%s)", yl, what),
                 caption = .wrap_caption(list(
                   c(paste0("How to read: each point is an individual once ", cov_txt, " have been removed from both variables: ",
                            "axes = residuals of a linear regression of ",
                            if (meth == "spearman") "the ranks of each variable on the ranks of the covariates" else "each variable on the covariates",
                            ". The correlation of these residuals is the partial correlation."),
                     fit$caption),
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
    gcol <- stats::setNames(rep_len(group_colours, 2), lev)
    rb_ci <- if (is.finite(eff$r_rb_low)) sprintf(" [%.2f, %.2f]", eff$r_rb_low, eff$r_rb_high) else ""
    eff_txt <- sprintf("r_rb = %.2f%s | mean difference %s \u2212 %s = %.2f [%.2f, %.2f] | max |r_s| for this split = %.2f",
                       eff$r_rb, rb_ci, lev[2], lev[1], eff$diff, eff$diff_low, eff$diff_high, eff$r_s_bound)
    return(ggplot2::ggplot(d, ggplot2::aes(x = g, y = v)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = g), width = 0.5, outlier.shape = NA, alpha = 0.25,
                            colour = "grey35", show.legend = FALSE) +
      ggplot2::geom_jitter(ggplot2::aes(colour = g), width = 0.12, height = 0, alpha = 0.7, size = 1.6,
                           show.legend = FALSE) +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3.4, fill = "white",
                            colour = "black", stroke = 0.9) +
      ggplot2::scale_fill_manual(values = gcol) + ggplot2::scale_colour_manual(values = gcol) +
      ggplot2::scale_x_discrete(labels = sprintf("%s\n(n = %d)", lev, c(eff$n0, eff$n1))) +
      .labs_safe(title = title %||% sprintf("%s by %s", vl, bl),
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
  measure <- if (meth == "pearson") "r measures how close the points are to a straight line"
             else sprintf("%s measures a monotonic (not necessarily linear) trend, computed on the ranks", sym)
  ggplot2::ggplot(d, ggplot2::aes(x = ex, y = ey)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8, colour = point_colour) +
    fit$layer +
    .labs_safe(title = title %||% sprintf("%s vs %s", yl, xl),
               subtitle = paste0(stat_txt, "\n", eng_txt), x = xl, y = yl,
               caption = .wrap_caption(list(
                 c(paste0("How to read: each point is an individual; ", measure, "."), fit$caption),
                 .cor_formula_lines(meth)), 110)) +
    thm
}

# Internal: fitted curve of plot_cor_pair() and the sentence that describes it
.smooth_layer <- function(smooth, smooth_args, colour, se) {
  if (identical(smooth, "none") || is.null(smooth)) return(list(layer = NULL, caption = NULL))
  name <- if (is.character(smooth)) smooth else "custom"
  args <- utils::modifyList(list(method = smooth, se = se, colour = colour, fill = colour, alpha = 0.15,
                                 linewidth = 1), smooth_args)
  if (is.null(args$formula) && name != "gam") args$formula <- y ~ x
  band <- if (se) " with its 95% confidence band" else ""
  txt <- switch(name,
                lm    = "Line = least-squares fit",
                loess = "Curve = LOESS smoother (local regression)",
                gam   = "Curve = GAM smoother (mgcv)",
                glm   = "Curve = generalised linear model fit",
                "Curve = fitted smoother")
  list(layer = do.call(ggplot2::geom_smooth, args),
       caption = paste0(txt, band, "; it shows the shape of the relation, the coefficient is in the subtitle."))
}
