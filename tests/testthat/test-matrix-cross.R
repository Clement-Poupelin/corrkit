d <- simulate_cor_data(150)
vars <- c("CRP", "IL6", "BDNF", "NfL", "BIN_CRP")

test_that("cor_test_matrix: symmetric, adjusted p >= raw p, engines agree", {
  cm <- cor_test_matrix(d, vars, engine = "stats", p_adjust = "BH")
  expect_s3_class(cm, "cor_matrix")
  expect_equal(cor_test_matrix(d, vars, engine = "stats")$p_adjust, "none")
  expect_equal(cm$r, t(cm$r))
  expect_equal(nrow(cm$long), choose(length(vars), 2))
  expect_true(all(cm$long$p_adj >= cm$long$p - 1e-12))
  expect_equal(unname(cm$types["BIN_CRP"]), "binary")
  skip_if_not_installed("psych")
  cp <- cor_test_matrix(d, vars, engine = "psych")
  expect_equal(cp$r, cm$r, tolerance = 1e-8)
  expect_equal(cp$long$p, cm$long$p, tolerance = 1e-6)
})

test_that("cor_test_matrix with the correlation engine", {
  skip_if_not_installed("correlation")
  cs <- cor_test_matrix(d, vars, engine = "stats")
  cc <- cor_test_matrix(d, vars, engine = "correlation")
  expect_equal(cc$r, cs$r, tolerance = 1e-8)
  ck <- cor_test_matrix(d, vars, method = "kendall")
  expect_equal(ck$method, "kendall")
})

test_that("cor_cross: adjustment within the chosen family", {
  cc <- cor_cross(d, c("CRP", "IL6", "BIN_CRP"), c("depression", "fatigue"), p_adjust = "BH", adjust_by = "y")
  expect_s3_class(cc, "cor_cross")
  tab <- cc$table
  for (s in unique(tab$Y)) {
    sub <- tab[tab$Y == s, ]
    expect_equal(sub$p_adj, stats::p.adjust(sub$p, "BH"))
  }
  none <- cor_cross(d, c("CRP", "IL6"), "fatigue", p_adjust = "BH", adjust_by = "none")
  expect_equal(none$table$p_adj, none$table$p)
  raw <- cor_cross(d, c("CRP", "IL6"), "fatigue")          # default: raw p-values
  expect_equal(raw$table$p_adj, raw$table$p)
  expect_equal(raw$meta$p_adjust, "none")
})

test_that("binary effects: r_s = r_rb * sqrt(3pq) without ties", {
  set.seed(1)
  dd <- data.frame(g = factor(rep(c("a", "b"), c(30, 70))), y = stats::rnorm(100))
  e <- cor_binary_effects(dd, "g", "y")
  expect_equal(e$r_s, e$r_rb * sqrt(3 * e$prop1 * (1 - e$prop1)), tolerance = 0.02)
  expect_true(abs(e$r_s) <= e$r_s_bound)
})
