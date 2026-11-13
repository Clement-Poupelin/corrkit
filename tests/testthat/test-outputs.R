d <- simulate_cor_data(120)
lab <- attr(d, "labels")

test_that("plots are ggplot objects with a caption", {
  cm <- cor_test_matrix(d, c("CRP", "IL6", "NfL", "BIN_CRP", "BIN_IL6"))
  p1 <- plot_cor_heatmap(cm, labels = lab)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(plot_cor_heatmap(cm, separator = "line"), "ggplot")
  res <- plot_cor_cross(d, c("CRP", "BIN_CRP"), c("depression", "fatigue"), covariates = "age", labels = lab)
  expect_s3_class(res$plot, "ggplot")
  expect_type(res$note, "character")
  expect_s3_class(plot_cor_pair(d, "IL6", "fatigue"), "ggplot")
  expect_s3_class(plot_cor_pair(d, "BIN_CRP", "depression"), "ggplot")
  expect_s3_class(plot_cor_pair(d, "NfL", "cognition", covariates = "age", meth = "pearson"), "ggplot")
})

test_that("tables are gt tables", {
  cm <- cor_test_matrix(d, c("CRP", "IL6", "BIN_CRP"))
  cc <- cor_cross(d, c("CRP", "IL6", "BIN_CRP"), c("depression", "fatigue"), covariates = "age", labels = lab)
  expect_s3_class(tbl_cor_matrix(cm, labels = lab), "gt_tbl")
  expect_s3_class(tbl_cor_long(cm), "gt_tbl")
  expect_s3_class(tbl_cor_cross(cc), "gt_tbl")
  expect_s3_class(tbl_cor_long(cc), "gt_tbl")
  expect_s3_class(tbl_cor_long(cor_cross(d, "CRP", c("depression", "fatigue"), p_adjust = "BH", adjust_by = "x")), "gt_tbl")
  res <- plot_cor_cross(cc)
  expect_s3_class(tbl_cor_cross(res), "gt_tbl")
  eff <- cor_binary_effects(d, "BIN_CRP", c("depression", "fatigue"))
  expect_s3_class(tbl_cor_binary(eff, labels = lab), "gt_tbl")
})

test_that("binary notes exist in both languages", {
  for (m in c("spearman", "pearson", "kendall")) for (l in c("en", "fr")) {
    expect_match(cor_binary_note(m, l), "doi:10.2466")
    expect_false(grepl("{MARK}", cor_binary_note(m, l, short = TRUE), fixed = TRUE))
  }
})

test_that("plots can be drawn on a pdf device (Latin-1 text only)", {
  cm <- cor_test_matrix(d, c("CRP", "IL6", "BIN_CRP"), method = "kendall")
  plots <- list(plot_cor_heatmap(cm),
                plot_cor_cross(d, c("CRP", "BIN_CRP"), "fatigue", covariates = "age")$plot,
                plot_cor_pair(d, "BIN_CRP", "depression"),
                plot_cor_pair(d, "NfL", "cognition", covariates = "age", meth = "spearman"))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in plots) expect_no_error(print(p))
})
