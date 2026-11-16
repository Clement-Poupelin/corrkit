d <- simulate_cor_data(120)
lab <- attr(d, "labels")
cm <- cor_compute(d, c("CRP", "IL6", "NfL", "BIN_CRP", "BIN_IL6"), labels = lab)
cc <- cor_compute(d, c("CRP", "IL6", "BIN_CRP"), c("depression", "fatigue"), covariates = "age",
                  p_adjust = "BH", labels = lab)

test_that("plot() returns a ggplot for every type and mode", {
  for (res in list(cm, cc)) {
    expect_s3_class(plot(res), "ggplot")
    expect_s3_class(plot(res, type = "forest"), "ggplot")
    expect_s3_class(suppressMessages(plot(res, type = "scatter")), "ggplot")
    expect_s3_class(autoplot(res), "ggplot")
  }
  expect_s3_class(plot(cc, type = "scatter", pair = c("IL6", "fatigue")), "ggplot")
  expect_s3_class(plot_cor_heatmap(cm, separator = "line", order = "original"), "ggplot")
  expect_s3_class(plot_cor_pair(d, "BIN_CRP", "depression"), "ggplot")
  expect_s3_class(plot_cor_pair(d, "NfL", "cognition", covariates = "age", meth = "pearson"), "ggplot")
})

test_that("tbl_cor() returns gt tables for every layout", {
  expect_s3_class(tbl_cor(cm), "gt_tbl")
  expect_s3_class(tbl_cor(cm, layout = "long", sort_by = "abs"), "gt_tbl")
  expect_s3_class(tbl_cor(cm, layout = "wide"), "gt_tbl")
  expect_s3_class(tbl_cor(cc), "gt_tbl")
  expect_s3_class(tbl_cor(cc, layout = "long"), "gt_tbl")
  expect_error(tbl_cor(cc, layout = "matrix"), "matrix analysis")
  expect_s3_class(tbl_cor(cor_compute(d, "CRP", c("depression", "fatigue"), p_adjust = "BH", adjust_by = "x"),
                          layout = "long"), "gt_tbl")
  expect_s3_class(tbl_cor_binary(cor_binary_effects(cc)), "gt_tbl")
})

test_that("print() shows the result", {
  expect_output(print(cm), "cor_result")
  expect_output(print(cc), "partial")
})

test_that("binary notes exist in both languages", {
  for (m in c("spearman", "pearson", "kendall")) for (l in c("en", "fr")) {
    expect_match(cor_binary_note(m, l), "doi:10.2466")
    expect_false(grepl("{MARK}", cor_binary_note(m, l, short = TRUE), fixed = TRUE))
  }
})

test_that("plots can be drawn on a pdf device (Latin-1 text only)", {
  ck <- cor_compute(d, c("CRP", "IL6", "BIN_CRP"), method = "kendall")
  plots <- list(plot(ck), plot(ck, type = "forest"), plot(cc), plot(cc, type = "forest"),
                plot_cor_pair(d, "BIN_CRP", "depression"),
                plot_cor_pair(d, "NfL", "cognition", covariates = "age", meth = "spearman"))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in plots) expect_no_error(print(p))
})

test_that("plain text mode (no markdown) also draws on a pdf device", {
  old <- options(corrkit.markdown = FALSE)
  on.exit(options(old), add = TRUE)
  p1 <- plot(cc, type = "forest")
  p2 <- plot_cor_pair(d, "BIN_CRP", "depression")
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(print(p1))
  expect_no_error(print(p2))
})

test_that("the notation is converted to subscripts, function names are left alone", {
  out <- corrkit:::.md_math("r_s = r_rb x sqrt(3pq); r(e_x.z, e_y.z); p < 0.05 * ; correlation::cor_to_ci()")
  expect_match(out, "<i>r</i><sub>s</sub>", fixed = TRUE)
  expect_match(out, "<i>r</i><sub>rb</sub>", fixed = TRUE)
  expect_match(out, "<i>e</i><sub>x.z</sub>", fixed = TRUE)
  expect_match(out, "&lt;", fixed = TRUE)
  expect_match(out, "\\*", fixed = TRUE)
  expect_match(out, "cor\\_to\\_ci", fixed = TRUE)
  expect_false(grepl("<sub>to", out, fixed = TRUE))
})

test_that("plot_cor_pair(): choice of the fitted curve", {
  for (s in list("lm", "loess", "glm", "none")) expect_s3_class(plot_cor_pair(d, "IL6", "fatigue", smooth = s), "ggplot")
  expect_s3_class(plot_cor_pair(d, "IL6", "fatigue", smooth = "loess", smooth_args = list(span = 0.5)), "ggplot")
  expect_s3_class(plot_cor_pair(d, "IL6", "fatigue", smooth = "lm", smooth_args = list(formula = y ~ poly(x, 2))), "ggplot")
  expect_s3_class(plot(cc, type = "scatter", pair = c("IL6", "fatigue"), smooth = "loess", se = FALSE), "ggplot")
  p <- plot_cor_pair(d, "IL6", "fatigue")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomSmooth"), logical(1))))
  p0 <- plot_cor_pair(d, "IL6", "fatigue", smooth = "none")
  expect_false(any(vapply(p0$layers, function(l) inherits(l$geom, "GeomSmooth"), logical(1))))
})
