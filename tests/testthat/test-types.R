test_that("variable types are detected", {
  expect_equal(cor_var_type(c(1.5, 2, 3.2)), "continuous")
  expect_equal(cor_var_type(c(0, 1, 1, NA)), "binary")
  expect_equal(cor_var_type(factor(c("low", "high"))), "binary")
  expect_equal(cor_var_type(factor(c("a", "b", "c"))), "nominal")
  expect_equal(cor_var_type(factor(c("a", "b", "c"), ordered = TRUE)), "ordinal")
})

test_that("binary coding puts the upper level at 1", {
  f <- factor(c("low", "high", "low"), levels = c("low", "high"))
  expect_equal(as_cor_numeric(f), c(0, 1, 0))
  expect_equal(as_cor_numeric(c(TRUE, FALSE)), c(1, 0))
  expect_equal(as_cor_numeric(c(1, 2, 2, NA)), c(0, 1, 1, NA))
  expect_warning(out <- as_cor_numeric(factor(c("a", "b", "c")), "v"), "unordered")
  expect_true(all(is.na(out)))
})

test_that("stars follow the thresholds", {
  expect_equal(p_to_stars(c(0.0001, 0.005, 0.03, 0.08, 0.5, NA), ns = "", trend = TRUE),
               c("***", "**", "*", ".", "", ""))
})
