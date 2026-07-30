.equa_methods <- c("kernel_end", "po", "kernel_ex")

.equa_method_label <- function(method) {
  switch(
    method,
    po = "Po-EQuA",
    kernel_end = "End-EQuA",
    kernel_ex = "Ex-EQuA",
    method
  )
}

.match_equa_method <- function(method) {
  match.arg(method, .equa_methods)
}

.validate_numeric_vector <- function(x, name, expected_length = NULL) {
  if (!is.numeric(x) || is.matrix(x) || is.data.frame(x) || is.complex(x)) {
    stop(sprintf("`%s` must be a numeric vector.", name), call. = FALSE)
  }
  if (!is.null(expected_length) && length(x) != expected_length) {
    stop(
      sprintf("`%s` must have length %d.", name, expected_length),
      call. = FALSE
    )
  }
  as.numeric(x)
}

.validate_probs <- function(probs) {
  probs <- .validate_numeric_vector(probs, "probs")
  if (!length(probs)) {
    stop("`probs` must contain at least one percentile level.", call. = FALSE)
  }
  if (any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("Every value in `probs` must be finite and strictly between 0 and 1.", call. = FALSE)
  }
  if (is.unsorted(probs, strictly = TRUE)) {
    stop("`probs` must be strictly increasing with no duplicates.", call. = FALSE)
  }
  probs
}

.prepare_covariates <- function(covariates, n, name = "covariates") {
  if (is.null(covariates)) {
    return(NULL)
  }

  if (is.data.frame(covariates)) {
    if (!ncol(covariates)) {
      stop(sprintf("`%s` must contain at least one column.", name), call. = FALSE)
    }
    numeric_columns <- vapply(covariates, is.numeric, logical(1))
    if (!all(numeric_columns)) {
      stop(sprintf("Every column in `%s` must be numeric.", name), call. = FALSE)
    }
    out <- as.matrix(covariates)
  } else if (is.matrix(covariates)) {
    if (!is.numeric(covariates) || is.complex(covariates)) {
      stop(sprintf("`%s` must be a numeric matrix or data frame.", name), call. = FALSE)
    }
    out <- covariates
  } else if (is.numeric(covariates) && !is.complex(covariates)) {
    out <- matrix(as.numeric(covariates), ncol = 1L)
  } else {
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  }

  if (nrow(out) != n) {
    stop(sprintf("`%s` must have %d rows.", name, n), call. = FALSE)
  }
  if (!ncol(out)) {
    stop(sprintf("`%s` must contain at least one column.", name), call. = FALSE)
  }

  column_names <- colnames(out)
  if (is.null(column_names)) {
    column_names <- paste0("z", seq_len(ncol(out)))
  }
  if (anyNA(column_names) || any(column_names == "") || anyDuplicated(column_names)) {
    stop(sprintf("`%s` must have unique, non-empty column names.", name), call. = FALSE)
  }
  colnames(out) <- column_names
  storage.mode(out) <- "double"
  out
}

.prepare_calibration_inputs <- function(
    chronological_age,
    estimated_age,
    covariates,
    method,
    na_action,
    warn = TRUE) {
  chronological_age <- .validate_numeric_vector(
    chronological_age,
    "chronological_age"
  )
  estimated_age <- .validate_numeric_vector(
    estimated_age,
    "estimated_age",
    expected_length = length(chronological_age)
  )
  n <- length(chronological_age)
  if (!n) {
    stop("Calibration inputs must not be empty.", call. = FALSE)
  }

  covariates <- .prepare_covariates(covariates, n)
  if (method == "kernel_ex" && is.null(covariates)) {
    stop("`covariates` are required for `method = \"kernel_ex\"`.", call. = FALSE)
  }
  if (method != "kernel_ex" && !is.null(covariates)) {
    stop(
      sprintf("`covariates` must be NULL for `method = \"%s\"`.", method),
      call. = FALSE
    )
  }

  complete <- is.finite(chronological_age) & is.finite(estimated_age)
  if (!is.null(covariates)) {
    complete <- complete & apply(is.finite(covariates), 1L, all)
  }

  n_omitted <- sum(!complete)
  if (n_omitted && na_action == "fail") {
    stop(
      sprintf("Calibration inputs contain %d incomplete or non-finite rows.", n_omitted),
      call. = FALSE
    )
  }
  if (n_omitted && isTRUE(warn)) {
    warning(
      sprintf("Omitted %d incomplete or non-finite calibration rows.", n_omitted),
      call. = FALSE
    )
  }
  if (sum(complete) < 3L) {
    stop("At least three complete calibration rows are required.", call. = FALSE)
  }

  list(
    n_input = n,
    complete = complete,
    n_omitted = n_omitted,
    chronological_age = chronological_age[complete],
    estimated_age = estimated_age[complete],
    covariates = if (is.null(covariates)) NULL else covariates[complete, , drop = FALSE]
  )
}

.validate_support_settings <- function(
    support_action,
    min_effective_n,
    nearest_k) {
  support_action <- match.arg(
    support_action,
    c("warn", "na", "error", "nearest")
  )
  min_effective_n <- .validate_numeric_vector(
    min_effective_n,
    "min_effective_n"
  )
  if (length(min_effective_n) != 1L ||
      !is.finite(min_effective_n) ||
      min_effective_n <= 0) {
    stop("`min_effective_n` must be one positive finite number.", call. = FALSE)
  }
  nearest_k <- .validate_numeric_vector(nearest_k, "nearest_k")
  if (length(nearest_k) != 1L ||
      !is.finite(nearest_k) ||
      nearest_k < 1 ||
      nearest_k != floor(nearest_k) ||
      nearest_k > .Machine$integer.max) {
    stop("`nearest_k` must be one positive integer.", call. = FALSE)
  }

  list(
    support_action = support_action,
    min_effective_n = as.numeric(min_effective_n),
    nearest_k = as.integer(nearest_k)
  )
}

.profile_column_names <- function(prefix, probs) {
  pct <- trimws(formatC(
    100 * probs,
    format = "fg",
    digits = 15,
    drop0trailing = TRUE
  ))
  pct <- gsub("\\.", "_", pct)
  pct <- gsub("-", "m", pct, fixed = TRUE)
  make.unique(paste0(prefix, "_p", pct), sep = "_")
}

.median_or_na <- function(x) {
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  stats::median(x[is.finite(x)])
}

.mean_or_na <- function(x) {
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  mean(x[is.finite(x)])
}
