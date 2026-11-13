# ====================================================================
# types.R -- variable types and numeric coding for a correlation
# ====================================================================

#' Type of a variable for a correlation
#'
#' @param x a vector (numeric, logical, factor or character)
#' @return `"binary"` (exactly 2 distinct non-missing values), `"ordinal"`
#'   (ordered factor with > 2 levels), `"nominal"` (unordered factor or
#'   character with > 2 levels) or `"continuous"`.
#' @export
#' @examples
#' cor_var_type(c(1.2, 3.4, 5.6))
#' cor_var_type(factor(c("low", "high", "low")))
#' cor_var_type(factor(c("a", "b", "c")))
cor_var_type <- function(x) {
  if (is.factor(x) || is.character(x)) {
    x <- droplevels(as.factor(x))
    if (nlevels(x) > 2) return(if (is.ordered(x)) "ordinal" else "nominal")
  }
  if (length(unique(stats::na.omit(x))) == 2) "binary" else "continuous"
}

#' Numeric coding of a variable for a correlation
#'
#' * logical -> 0/1;
#' * binary factor or character -> 0/1, the upper level coded 1 (so r > 0
#'   means "higher in the upper level"); binary numeric -> 0/1 (largest value = 1);
#' * ordered factor with > 2 levels -> 0, 1, 2, ... (Pearson then assumes
#'   equal spacing between levels: warning);
#' * unordered factor with > 2 levels -> NA with a warning (no meaningful correlation).
#'
#' Recoding a binary variable to 0/1 is a positive linear transformation: it
#' changes neither r, nor r_s, nor tau, nor their p-values.
#' @param x vector to code
#' @param name name used in the warnings
#' @param meth correlation method (`"pearson"` triggers the equal-spacing warning)
#' @return numeric vector
#' @export
#' @examples
#' as_cor_numeric(factor(c("low", "high", "low"), levels = c("low", "high")))
as_cor_numeric <- function(x, name = "x", meth = "spearman") {
  if (is.logical(x)) return(as.numeric(x))
  if (is.factor(x) || is.character(x)) {
    x <- droplevels(as.factor(x))
    k <- nlevels(x)
    if (k > 2 && !is.ordered(x)) {
      warning(sprintf("%s: unordered factor with %d levels, not usable in a correlation (set to NA).", name, k),
              call. = FALSE)
      return(rep(NA_real_, length(x)))
    }
    if (k > 2 && meth == "pearson")
      warning(sprintf("%s: ordered factor coded 0..%d; Pearson assumes equal spacing between levels.", name, k - 1),
              call. = FALSE)
    return(as.numeric(x) - 1)
  }
  x <- as.numeric(x)
  u <- sort(unique(stats::na.omit(x)))
  if (length(u) == 2) x <- as.numeric(x == u[2])
  x
}

#' Add a marker (default `" [0/1]"`) to the labels of binary variables
#'
#' @param lbl labels
#' @param types types from [cor_var_type()], same length as `lbl`
#' @param mark marker appended to the binary labels (explained in captions and notes)
#' @export
mark_binary <- function(lbl, types, mark = " [0/1]") ifelse(types == "binary", paste0(lbl, mark), lbl)
