d <- simulate_cor_data(150)
vars <- c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP")

test_that("matrix mode: one row per pair, symmetric matrix, raw p by default", {
  cm <- cor_compute(d, vars, engine = "stats")
  expect_s3_class(cm, "cor_result")
  expect_equal(cm$meta$mode, "matrix")
  expect_equal(nrow(cm$table), choose(length(vars), 2))
  R <- as.matrix(cm)
  expect_equal(R, t(R))
  expect_equal(unname(diag(R)), rep(1, length(vars)))
  expect_equal(cm$meta$p_adjust, "none")
  expect_equal(cm$table$p_adj, cm$table$p)
  expect_equal(unname(cm$meta$x_types["BIN_CRP"]), "binary")
})

test_that("matrix mode: engines agree", {
  cs <- cor_compute(d, vars, engine = "stats")
  skip_if_not_installed("psych")
  cp <- cor_compute(d, vars, engine = "psych")
  expect_equal(as.matrix(cp), as.matrix(cs), tolerance = 1e-8)
  expect_equal(cp$table$p, cs$table$p, tolerance = 1e-6)
})

test_that("matrix mode: correlation engine, Kendall and robust methods", {
  skip_if_not_installed("correlation")
  cs <- cor_compute(d, vars, engine = "stats")
  cc <- cor_compute(d, vars, engine = "correlation")
  expect_equal(as.matrix(cc), as.matrix(cs), tolerance = 1e-8)
  expect_equal(cor_compute(d, vars, method = "kendall")$meta$meth, "kendall")
  expect_s3_class(cor_compute(d, vars[1:3], method = "biweight"), "cor_result")
})

test_that("matrix and pair-by-pair computations give the same result", {
  a <- cor_compute(d, vars, engine = "stats")
  b <- cor_compute(d, vars[1:2], vars[3:5])
  ab <- merge(as.data.frame(a)[, c("X", "Y", "r", "p")], as.data.frame(b)[, c("X", "Y", "r", "p")],
              by = c("X", "Y"))
  expect_true(nrow(ab) > 0)
  expect_equal(ab$r.x, ab$r.y, tolerance = 1e-8)
  expect_equal(ab$p.x, ab$p.y, tolerance = 1e-6)
})

test_that("cross mode: adjustment within the chosen family", {
  cc <- cor_compute(d, c("CRP", "IL6", "BIN_CRP"), c("depression", "fatigue"), p_adjust = "BH", adjust_by = "y")
  expect_equal(cc$meta$mode, "cross")
  expect_equal(dim(as.matrix(cc)), c(3L, 2L))
  for (s in unique(cc$table$Y)) {
    sub <- cc$table[cc$table$Y == s, ]
    expect_equal(sub$p_adj, stats::p.adjust(sub$p, "BH"))
  }
  all <- cor_compute(d, c("CRP", "IL6"), c("depression", "fatigue"), p_adjust = "holm", adjust_by = "all")
  expect_equal(all$table$p_adj, stats::p.adjust(all$table$p, "holm"))
  raw <- cor_compute(d, c("CRP", "IL6"), "fatigue")
  expect_equal(raw$table$p_adj, raw$table$p)
})

test_that("partial correlations: matrix and cross modes, equal to partial_cor_test()", {
  pm <- cor_compute(d, c("NfL", "cognition", "CRP"), covariates = c("age", "sex"), method = "pearson")
  expect_true(pm$meta$partial)
  expect_equal(pm$meta$k, 2L)
  one <- partial_cor_test(d$NfL, d$cognition, d[, c("age", "sex")], meth = "pearson")
  row <- pm$table[pm$table$X == "NfL" & pm$table$Y == "cognition", ]
  expect_equal(row$r, one$estimate, tolerance = 1e-8)
  expect_equal(row$p, one$p.value, tolerance = 1e-8)
})

test_that("single pair", {
  one <- cor_compute(d, "IL6", "fatigue")
  expect_equal(nrow(one$table), 1L)
})

test_that("summary and as.data.frame", {
  s <- summary(cor_compute(d, vars, p_adjust = "BH"))
  expect_s3_class(s, "summary.cor_result")
  expect_equal(nrow(s$variables), length(vars))
  expect_output(print(s), "Correlation analysis")
  expect_s3_class(as.data.frame(cor_compute(d, vars)), "data.frame")
})

test_that("binary effects: r_s = r_rb * sqrt(3pq) without ties; works from a cor_result", {
  set.seed(1)
  dd <- data.frame(g = factor(rep(c("a", "b"), c(30, 70))), y = stats::rnorm(100))
  e <- cor_binary_effects(dd, "g", "y")
  expect_equal(e$r_s, e$r_rb * sqrt(3 * e$prop1 * (1 - e$prop1)), tolerance = 0.02)
  expect_true(abs(e$r_s) <= e$r_s_bound)
  e2 <- cor_binary_effects(cor_compute(d, c("BIN_CRP", "CRP"), c("depression", "fatigue")))
  expect_equal(sort(unique(e2$X)), "BIN_CRP")
  expect_equal(nrow(e2), 2L)
})
