.weighted_quantile_left <- function(x, weights, probs) {
  ok <- is.finite(x) & is.finite(weights) & weights > 0
  x <- x[ok]
  weights <- weights[ok]

  if (!length(x)) {
    return(rep(NA_real_, length(probs)))
  }

  order_index <- order(x)
  x <- x[order_index]
  weights <- weights[order_index]

  # Scaling by the maximum weight prevents overflow without changing the
  # weighted empirical distribution.
  weights <- weights / max(weights)
  total_weight <- sum(weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    return(rep(NA_real_, length(probs)))
  }

  cumulative_weight <- cumsum(weights / total_weight)
  vapply(
    probs,
    function(prob) {
      index <- which(cumulative_weight >= prob)[1L]
      if (is.na(index)) x[[length(x)]] else x[[index]]
    },
    numeric(1)
  )
}
