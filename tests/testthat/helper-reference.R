reference_weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x) || sum(w) <= 0) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  cumulative <- cumsum(w)
  vapply(
    probs,
    function(prob) x[which(cumulative >= prob)[1L]],
    numeric(1)
  )
}

reference_bandwidth <- function(x, dimension) {
  x <- x[is.finite(x)]
  max(
    1.06 * stats::sd(x) * length(x)^(-1 / (4 + dimension)),
    0.25 * stats::IQR(x),
    0.05 * diff(range(x)),
    1e-6,
    na.rm = TRUE
  )
}

reference_age_folds <- function(age, folds, seed, n_strata = 10L) {
  n <- length(age)
  n_strata_use <- min(n_strata, max(1L, floor(n / folds)))
  age_order <- order(age, seq_along(age))
  stratum <- integer(n)
  stratum[age_order] <-
    ((seq_len(n) - 1L) * n_strata_use) %/% n + 1L

  set.seed(seed)
  assignment <- integer(n)
  for (stratum_id in seq_len(n_strata_use)) {
    index <- age_order[stratum[age_order] == stratum_id]
    shuffled <- sample(index, length(index))
    assignment[shuffled] <- rep(seq_len(folds), length.out = length(index))
  }
  assignment
}
