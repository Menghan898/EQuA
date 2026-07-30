#' Extract EQuA prediction diagnostics
#'
#' @param x An `equa_profile` or `equa_crossfit` object.
#' @param by_variable Whether to return marginal range diagnostics for every
#'   conditioning variable.
#'
#' @return A data frame, or a list containing summary and by-variable
#'   diagnostics.
#'
#' @examples
#' calibration <- subset(
#'   simulated_equa_data,
#'   sample_role == "calibration"
#' )
#' fit <- fit_equa(
#'   calibration$chronological_age,
#'   calibration$estimated_age,
#'   method = "kernel_end"
#' )
#' profile <- predict(
#'   fit,
#'   new_estimated_age = calibration$estimated_age[1:3]
#' )
#' equa_diagnostics(profile)
#'
#' @export
equa_diagnostics <- function(x, by_variable = FALSE) {
  if (!inherits(x, "equa_profile")) {
    stop("`x` must inherit from `equa_profile`.", call. = FALSE)
  }
  if (!is.logical(by_variable) || length(by_variable) != 1L || is.na(by_variable)) {
    stop("`by_variable` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!by_variable) {
    return(x$diagnostics)
  }
  list(
    summary = x$diagnostics,
    outside_range_by_variable = x$outside_range_by_variable
  )
}

#' @export
print.equa_fit <- function(x, ...) {
  cat("<equa_fit>\n")
  cat("  Method:               ", .equa_method_label(x$method), "\n", sep = "")
  cat("  Calibration rows:     ", x$n_calibration, "\n", sep = "")
  cat("  Omitted rows:         ", x$n_omitted, "\n", sep = "")
  cat("  Percentile levels:    ", length(x$probs), "\n", sep = "")
  cat("  Conditioning dimensions: ", length(x$conditioning_names), "\n", sep = "")
  cat("  Support action:       ", x$support_action, "\n", sep = "")
  cat("  Privacy: fitted object retains calibration-level values\n")
  invisible(x)
}

#' @export
summary.equa_fit <- function(object, ...) {
  out <- list(
    method = object$method,
    n_input = object$n_input,
    n_calibration = object$n_calibration,
    n_omitted = object$n_omitted,
    probs = object$probs,
    conditioning_names = object$conditioning_names,
    bandwidth = object$bandwidth,
    bandwidth_selected = object$bandwidth_selected,
    na_action = object$na_action,
    support_action = object$support_action,
    min_effective_n = object$min_effective_n,
    privacy_notice = object$privacy_notice
  )
  class(out) <- "summary_equa_fit"
  out
}

#' @export
print.summary_equa_fit <- function(x, ...) {
  cat("EQuA fit summary\n")
  cat("  Method:           ", .equa_method_label(x$method), "\n", sep = "")
  cat("  Calibration rows: ", x$n_calibration, "\n", sep = "")
  cat("  Omitted rows:     ", x$n_omitted, "\n", sep = "")
  cat("  Probabilities:    ", paste(x$probs, collapse = ", "), "\n", sep = "")
  if (length(x$bandwidth)) {
    cat(
      "  Bandwidths:       ",
      paste(
        paste0(names(x$bandwidth), "=", format(x$bandwidth, digits = 4)),
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }
  cat("  Privacy: fitted object retains calibration-level values\n")
  invisible(x)
}

#' @export
print.equa_profile <- function(x, ...) {
  cat("<", class(x)[[1L]], ">\n", sep = "")
  cat("  Method:            ", .equa_method_label(x$method), "\n", sep = "")
  cat("  Observations:      ", nrow(x$qea), "\n", sep = "")
  cat("  Percentile levels: ", ncol(x$qea), "\n", sep = "")
  cat(
    "  Low-support rows:  ",
    sum(x$diagnostics$low_support, na.rm = TRUE),
    "\n",
    sep = ""
  )
  cat(
    "  Incomplete rows:   ",
    sum(x$diagnostics$incomplete_conditioning, na.rm = TRUE),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
print.equa_crossfit <- function(x, ...) {
  NextMethod("print")
  cat("  Folds:             ", x$config$folds, "\n", sep = "")
  cat("  Repeats:           ", x$config$repeats, "\n", sep = "")
  cat("  Aggregation:       ", x$config$aggregate, "\n", sep = "")
  invisible(x)
}

#' @export
summary.equa_profile <- function(object, ...) {
  out <- list(
    class = class(object)[[1L]],
    method = object$method,
    n = nrow(object$qea),
    probs = object$probs,
    qea_available = rowSums(is.finite(object$qea)) > 0,
    qeaa_available = if (is.null(object$qeaa)) {
      rep(FALSE, nrow(object$qea))
    } else {
      rowSums(is.finite(object$qeaa)) > 0
    },
    n_low_support = sum(object$diagnostics$low_support, na.rm = TRUE),
    n_outside_training_range = sum(
      object$diagnostics$outside_training_range,
      na.rm = TRUE
    ),
    n_incomplete = sum(
      object$diagnostics$incomplete_conditioning,
      na.rm = TRUE
    )
  )
  if (inherits(object, "equa_crossfit")) {
    out$config <- object$config
  }
  class(out) <- "summary_equa_profile"
  out
}

#' @export
summary.equa_crossfit <- function(object, ...) {
  NextMethod("summary")
}

#' @export
print.summary_equa_profile <- function(x, ...) {
  cat("EQuA profile summary\n")
  cat("  Object class:       ", x$class, "\n", sep = "")
  cat("  Method:             ", .equa_method_label(x$method), "\n", sep = "")
  cat("  Observations:       ", x$n, "\n", sep = "")
  cat("  qEA available:      ", sum(x$qea_available), "\n", sep = "")
  cat("  qEAA available:     ", sum(x$qeaa_available), "\n", sep = "")
  cat("  Low-support rows:   ", x$n_low_support, "\n", sep = "")
  cat("  Outside-range rows: ", x$n_outside_training_range, "\n", sep = "")
  cat("  Incomplete rows:    ", x$n_incomplete, "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.equa_profile <- function(
    x,
    row.names = NULL,
    optional = FALSE,
    ...,
    format = c("wide", "long")) {
  format <- match.arg(format)
  n <- nrow(x$qea)
  base <- data.frame(
    observation = seq_len(n),
    estimated_age = x$estimated_age,
    stringsAsFactors = FALSE
  )
  if (!is.null(x$chronological_age)) {
    base$chronological_age <- x$chronological_age
  }

  if (format == "wide") {
    out <- cbind(base, as.data.frame(x$qea, optional = TRUE))
    if (!is.null(x$qeaa)) {
      out <- cbind(out, as.data.frame(x$qeaa, optional = TRUE))
    }
    return(cbind(out, x$diagnostics))
  }

  probability_index <- rep(seq_along(x$probs), each = n)
  observation_index <- rep(seq_len(n), times = length(x$probs))
  out <- base[observation_index, , drop = FALSE]
  out$tau <- x$probs[probability_index]
  out$qea <- as.vector(x$qea)
  out$qeaa <- if (is.null(x$qeaa)) NA_real_ else as.vector(x$qeaa)
  diagnostics <- x$diagnostics[observation_index, , drop = FALSE]
  rownames(out) <- NULL
  cbind(out, diagnostics)
}

#' @export
as.data.frame.equa_crossfit <- function(x, ...) {
  NextMethod("as.data.frame")
}

#' Plot EQuA profiles
#'
#' @param x An `equa_profile` or `equa_crossfit` object.
#' @param type Either `"qea"` or `"qeaa"`.
#' @param rows Integer row indices to display.
#' @param xlab,ylab Axis labels.
#' @param ... Additional arguments passed to [graphics::matplot()].
#'
#' @return The plotted object, invisibly.
#' @export
plot.equa_profile <- function(
    x,
    type = c("qea", "qeaa"),
    rows = seq_len(min(6L, nrow(x$qea))),
    xlab = "Percentile level",
    ylab = NULL,
    ...) {
  type <- match.arg(type)
  values <- if (type == "qea") x$qea else x$qeaa
  if (is.null(values)) {
    stop("qEAA values are unavailable because chronological age was not supplied.", call. = FALSE)
  }
  if (!is.numeric(rows) ||
      !length(rows) ||
      any(!is.finite(rows)) ||
      any(rows != floor(rows)) ||
      any(rows < 1 | rows > nrow(values))) {
    stop("`rows` must contain valid integer observation indices.", call. = FALSE)
  }
  rows <- as.integer(rows)
  selected <- values[rows, , drop = FALSE]
  if (!any(is.finite(selected))) {
    stop("The selected profiles contain no finite values.", call. = FALSE)
  }
  if (is.null(ylab)) {
    ylab <- if (type == "qea") "Quantile age" else "Quantile age acceleration"
  }

  graphics::matplot(
    x$probs,
    t(selected),
    type = "l",
    lty = 1,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  invisible(x)
}

#' @export
plot.equa_crossfit <- function(x, ...) {
  NextMethod("plot")
}
