# ====================================================================
# simulate.R -- example data used in the examples, tests and vignettes
# ====================================================================

#' Simulated biomarkers, clinical scores and covariates
#'
#' Data generated on the fly (no stored data set), with known structure:
#' * a latent inflammation L drives `CRP` (log-normal, skewed), `IL6`,
#'   `BDNF` (negatively) and the `depression` score;
#' * `IL6` drives `fatigue`;
#' * `age` drives both `NfL` and `cognition`: their raw correlation is mostly
#'   confounded by age and shrinks once adjusted for age;
#' * `BIN_CRP` = CRP above its median ("high" vs "low", 50/50 split) and
#'   `BIN_IL6` = IL6 above its 80th percentile (20/80 split): the same
#'   information as the continuous markers, cut into two groups;
#' * about `na_rate` missing values in the biomarkers and scores.
#'
#' The random-number state of the session is restored on exit.
#' @param n number of individuals
#' @param seed seed of the simulation
#' @param na_rate share of missing values in each biomarker and score
#' @return data.frame with an attribute `"labels"` (named vector variable -> label)
#' @export
#' @examples
#' d <- simulate_cor_data(200)
#' str(d)
#' attr(d, "labels")
simulate_cor_data <- function(n = 200, seed = 2026, na_rate = 0.03) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  rn  <- function(s = 1) stats::rnorm(n, 0, s)
  age <- round(stats::rnorm(n, 55, 12))
  sex <- factor(sample(c("F", "M"), n, replace = TRUE), levels = c("F", "M"))
  L   <- rn()
  a   <- (age - 55) / 12
  d <- data.frame(
    age = age, sex = sex,
    CRP  = round(exp(0.6 * L + 0.1 * a + rn(0.6)), 2),
    IL6  = round(0.6 * L + rn(0.8), 2),
    BDNF = round(-0.4 * L + 0.3 * (sex == "M") + rn(), 2),
    NfL  = round(0.7 * a + rn(0.7), 2)
  )
  d$depression <- round(12 + 2.0 * L + 0.4 * a + rn(3), 1)
  d$fatigue    <- round(8 + 1.5 * d$IL6 + rn(3), 1)
  d$cognition  <- round(26 - 1.2 * a - 0.3 * d$NfL + rn(1.5), 1)
  d$BIN_CRP <- factor(ifelse(d$CRP > stats::median(d$CRP), "high", "low"), levels = c("low", "high"))
  d$BIN_IL6 <- factor(ifelse(d$IL6 > stats::quantile(d$IL6, 0.8), "high", "low"), levels = c("low", "high"))
  if (na_rate > 0) for (v in c("CRP", "IL6", "BDNF", "NfL", "depression", "fatigue", "cognition")) {
    d[[v]][sample.int(n, max(0, round(na_rate * n)))] <- NA
  }
  # binary versions keep the cut-off of the complete marker; missing where the marker is missing
  d$BIN_CRP[is.na(d$CRP)] <- NA
  d$BIN_IL6[is.na(d$IL6)] <- NA
  attr(d, "labels") <- c(age = "Age", sex = "Sex", CRP = "CRP", IL6 = "IL-6", BDNF = "BDNF", NfL = "NfL",
                         depression = "Depression score", fatigue = "Fatigue score", cognition = "Cognition score",
                         BIN_CRP = "CRP high", BIN_IL6 = "IL-6 high")
  d
}
