d <- simulate_cor_data(150)

test_that("unadjusted test equals stats::cor.test", {
  ok <- stats::complete.cases(d$IL6, d$fatigue)
  ref <- suppressWarnings(stats::cor.test(d$IL6[ok], d$fatigue[ok], method = "spearman", exact = FALSE))
  res <- partial_cor_test(d$IL6, d$fatigue, meth = "spearman")
  expect_equal(res$estimate, unname(ref$estimate))
  expect_equal(res$p.value, ref$p.value)
  expect_equal(res$k, 0L)
  expect_true(res$ci_low < res$estimate && res$estimate < res$ci_high)
})

test_that("Pearson partial = slope test of lm (Frisch-Waugh-Lovell), both engines", {
  Z <- d[, c("age", "sex")]
  fit <- stats::lm(cognition ~ NfL + age + sex, data = d)
  p_lm <- summary(fit)$coefficients["NfL", "Pr(>|t|)"]
  base <- partial_cor_test(d$NfL, d$cognition, Z, meth = "pearson", engine = "base")
  expect_equal(base$p.value, p_lm, tolerance = 1e-8)
  expect_equal(base$k, 2L)
  skip_if_not_installed("ppcor")
  pp <- partial_cor_test(d$NfL, d$cognition, Z, meth = "pearson", engine = "ppcor")
  expect_equal(pp$estimate, base$estimate, tolerance = 1e-8)
  expect_equal(pp$p.value, base$p.value, tolerance = 1e-6)
})

test_that("Spearman partial: residual method agrees with ppcor", {
  skip_if_not_installed("ppcor")
  Z <- d[, "age", drop = FALSE]
  a <- partial_cor_test(d$CRP, d$depression, Z, meth = "spearman", engine = "ppcor")
  b <- partial_cor_test(d$CRP, d$depression, Z, meth = "spearman", engine = "base")
  expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
})

test_that("point-biserial test = Student's t-test", {
  g <- as_cor_numeric(d$BIN_CRP); y <- d$depression
  ok <- stats::complete.cases(g, y)
  res <- partial_cor_test(g, y, meth = "pearson")
  tt <- stats::t.test(y[ok & g == 1], y[ok & g == 0], var.equal = TRUE)
  expect_equal(res$p.value, tt$p.value, tolerance = 1e-10)
})
