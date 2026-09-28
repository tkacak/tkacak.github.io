threshold_table <- list(
  symmetric = list(
    `2` = 0,
    `3` = c(-0.83, 0.83),
    `4` = c(-1.25, 0.00, 1.25),
    `5` = c(-1.50, -0.50, 0.50, 1.50),
    `6` = c(-1.60, -0.83, 0.00, 0.83, 1.60),
    `7` = c(-1.79, -1.07, -0.36, 0.36, 1.07, 1.79)),
  mod_asym = list(
    `2` = 0.36,
    `3` = c(-0.50, 0.76),
    `4` = c(-0.31, 0.79, 1.66),
    `5` = c(-0.70, 0.39, 1.16, 2.05),
    `6` = c(-1.05, 0.08, 0.81, 1.44, 2.33),
    `7` = c(-1.43, -0.43, 0.38, 0.94, 1.44, 2.54)),
  ext_asym = list(
    `2` = 1.04,
    `3` = c(0.58, 1.13),
    `4` = c(0.28, 0.71, 1.23),
    `5` = c(0.05, 0.44, 0.84, 1.34),
    `6` = c(-0.13, 0.25, 0.61, 0.99, 1.48),
    `7` = c(-0.25, 0.13, 0.47, 0.81, 1.18, 1.64)),
  mod_asym_alt = list(
    `2` = -0.36,
    `3` = c(-0.76, 0.50),
    `4` = c(-1.66, -0.79, 0.31),
    `5` = c(-2.05, -1.16, -0.39, 0.70),
    `6` = c(-2.33, -1.44, -0.81, -0.08, 1.05),
    `7` = c(-2.54, -1.44, -0.94, -0.38, 0.43, 1.43)),
  ext_asym_alt = list(
    `2` = -1.04,
    `3` = c(-1.13, -0.58),
    `4` = c(-1.23, -0.71, -0.28),
    `5` = c(-1.34, -0.84, -0.44, -0.05),
    `6` = c(-1.48, -0.99, -0.61, -0.25, 0.13),
    `7` = c(-1.64, -1.18, -0.81, -0.47, -0.13, 0.25))
)

get_thresholds <- function(ncat, condition = "symmetric") {
  condition <- match.arg(condition, names(threshold_table))
  tau <- threshold_table[[condition]][[as.character(ncat)]]
  if (is.null(tau)) stop("ncat must be between 2 and 7, got ", ncat)
  tau
}

categorize_thresholds <- function(data, ncat, condition = "symmetric", standardize = FALSE) {
  X <- as.matrix(data)
  if (!is.numeric(X)) stop("data must be numeric")
  p <- ncol(X)
  ncat      <- rep_len(ncat, p)
  condition <- rep_len(condition, p)
  if (standardize) X <- scale(X)
  out <- vapply(seq_len(p), function(j)
    findInterval(X[, j], get_thresholds(ncat[j], condition[j])) + 1L,
    integer(nrow(X)))
  out <- matrix(out, nrow(X), p, dimnames = dimnames(as.matrix(data)))
  out
}
