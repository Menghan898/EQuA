.select_bandwidth_rule <- function(x, dimension) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L) {
    stop("At least three finite values are required to select a bandwidth.", call. = FALSE)
  }
  if (length(dimension) != 1L ||
      !is.finite(dimension) ||
      dimension < 1 ||
      dimension != floor(dimension)) {
    stop("`dimension` must be one positive integer.", call. = FALSE)
  }

  exponent <- 1 / (4 + dimension)
  bandwidth <- max(
    1.06 * stats::sd(x) * n^(-exponent),
    0.25 * stats::IQR(x),
    0.05 * diff(range(x)),
    1e-6,
    na.rm = TRUE
  )
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop("Could not determine a positive bandwidth.", call. = FALSE)
  }
  bandwidth
}

.conditioning_transform <- function(conditioning) {
  if (!ncol(conditioning)) {
    return(list(
      raw = conditioning,
      center = numeric(),
      scale = numeric(),
      scaled = conditioning,
      range = matrix(
        numeric(),
        nrow = 2L,
        ncol = 0L,
        dimnames = list(c("minimum", "maximum"), character())
      )
    ))
  }

  center <- colMeans(conditioning)
  scale <- apply(conditioning, 2L, stats::sd)
  invalid_scale <- !is.finite(scale) | scale <= 0
  if (any(invalid_scale)) {
    bad <- paste(colnames(conditioning)[invalid_scale], collapse = ", ")
    stop(
      sprintf("Conditioning variables must have positive variation: %s.", bad),
      call. = FALSE
    )
  }

  scaled <- sweep(conditioning, 2L, center, "-")
  scaled <- sweep(scaled, 2L, scale, "/")
  raw_range <- rbind(
    minimum = apply(conditioning, 2L, min),
    maximum = apply(conditioning, 2L, max)
  )

  list(
    raw = conditioning,
    center = center,
    scale = scale,
    scaled = scaled,
    range = raw_range
  )
}

.resolve_bandwidth <- function(bandwidth, transform) {
  dimension <- ncol(transform$scaled)
  if (!dimension) {
    if (!is.null(bandwidth)) {
      stop("`bandwidth` must be NULL for Po-EQuA.", call. = FALSE)
    }
    return(list(raw = numeric(), scaled = numeric(), selected = FALSE))
  }

  if (is.null(bandwidth)) {
    # Select bandwidths in the original variable units, exactly as in the
    # reproducibility programs. Convert them only after selection so the
    # internally standardized kernel evaluates the same raw-unit ratios.
    raw <- vapply(
      seq_len(dimension),
      function(j) .select_bandwidth_rule(transform$raw[, j], dimension),
      numeric(1)
    )
    scaled <- raw / transform$scale
    names(raw) <- colnames(transform$scaled)
    names(scaled) <- colnames(transform$scaled)
    return(list(raw = raw, scaled = scaled, selected = TRUE))
  }

  if (!is.numeric(bandwidth) || is.complex(bandwidth)) {
    stop("`bandwidth` must be a positive numeric vector.", call. = FALSE)
  }
  bandwidth <- as.numeric(bandwidth)
  if (length(bandwidth) != dimension) {
    stop(
      sprintf("`bandwidth` must have length %d.", dimension),
      call. = FALSE
    )
  }
  if (any(!is.finite(bandwidth)) || any(bandwidth <= 0)) {
    stop("Every bandwidth must be positive and finite.", call. = FALSE)
  }

  names(bandwidth) <- colnames(transform$scaled)
  scaled <- bandwidth / transform$scale
  names(scaled) <- names(bandwidth)
  list(raw = bandwidth, scaled = scaled, selected = FALSE)
}
