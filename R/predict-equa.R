.prepare_prediction_conditioning <- function(
    object,
    new_estimated_age,
    new_covariates) {
  new_estimated_age <- .validate_numeric_vector(
    new_estimated_age,
    "new_estimated_age"
  )
  n <- length(new_estimated_age)
  if (!n) {
    stop("Prediction inputs must not be empty.", call. = FALSE)
  }

  if (object$method == "kernel_ex") {
    new_covariates <- .prepare_covariates(
      new_covariates,
      n,
      name = "new_covariates"
    )
    expected_names <- object$conditioning_names[-1L]
    if (ncol(new_covariates) != length(expected_names)) {
      stop(
        sprintf(
          "`new_covariates` must contain %d columns.",
          length(expected_names)
        ),
        call. = FALSE
      )
    }
    if (!identical(colnames(new_covariates), expected_names)) {
      stop(
        paste(
          "`new_covariates` column names and order must match the",
          "calibration covariates."
        ),
        call. = FALSE
      )
    }
    raw <- cbind(estimated_age = new_estimated_age, new_covariates)
  } else {
    if (!is.null(new_covariates)) {
      stop(
        sprintf(
          "`new_covariates` must be NULL for method \"%s\".",
          object$method
        ),
        call. = FALSE
      )
    }
    raw <- if (object$method == "kernel_end") {
      matrix(
        new_estimated_age,
        ncol = 1L,
        dimnames = list(NULL, "estimated_age")
      )
    } else {
      matrix(
        numeric(),
        nrow = n,
        ncol = 0L,
        dimnames = list(NULL, character())
      )
    }
  }

  valid <- is.finite(new_estimated_age)
  if (ncol(raw)) {
    valid <- valid & apply(is.finite(raw), 1L, all)
  }

  scaled <- raw
  if (ncol(raw)) {
    scaled <- sweep(raw, 2L, object$conditioning_center, "-")
    scaled <- sweep(scaled, 2L, object$conditioning_scale, "/")
  }

  list(
    estimated_age = new_estimated_age,
    raw = raw,
    scaled = scaled,
    valid = valid
  )
}

.kernel_prediction_row <- function(object, new_scaled) {
  difference <- sweep(
    object$calibration_conditioning,
    2L,
    new_scaled,
    "-"
  )
  standardized_kernel_distance <- sweep(
    difference,
    2L,
    object$bandwidth_scaled,
    "/"
  )
  log_weights <- -0.5 * rowSums(standardized_kernel_distance^2)

  finite_log_weights <- is.finite(log_weights)
  if (!any(finite_log_weights)) {
    return(list(
      quantiles = rep(NA_real_, length(object$probs)),
      effective_sample_size = NA_real_,
      maximum_normalized_weight = NA_real_,
      normalization_failed = TRUE,
      nearest_standardized_distance = NA_real_,
      normalized_weights = NULL,
      euclidean_distances = rep(NA_real_, object$n_calibration)
    ))
  }

  stabilized <- exp(log_weights - max(log_weights[finite_log_weights]))
  stabilized[!finite_log_weights] <- 0
  total <- sum(stabilized)
  if (!is.finite(total) || total <= 0) {
    return(list(
      quantiles = rep(NA_real_, length(object$probs)),
      effective_sample_size = NA_real_,
      maximum_normalized_weight = NA_real_,
      normalization_failed = TRUE,
      nearest_standardized_distance = NA_real_,
      normalized_weights = NULL,
      euclidean_distances = rep(NA_real_, object$n_calibration)
    ))
  }

  normalized_weights <- stabilized / total
  euclidean_distances <- sqrt(rowSums(difference^2))
  list(
    quantiles = .weighted_quantile_left(
      object$calibration_ead,
      normalized_weights,
      object$probs
    ),
    effective_sample_size = 1 / sum(normalized_weights^2),
    maximum_normalized_weight = max(normalized_weights),
    normalization_failed = FALSE,
    nearest_standardized_distance = min(euclidean_distances),
    normalized_weights = normalized_weights,
    euclidean_distances = euclidean_distances
  )
}

.summarize_support_warning <- function(diagnostics) {
  n_low <- sum(diagnostics$low_support, na.rm = TRUE)
  n_failed <- sum(diagnostics$weight_normalization_failed, na.rm = TRUE)
  n_incomplete <- sum(diagnostics$incomplete_conditioning, na.rm = TRUE)
  parts <- character()
  if (n_low) {
    parts <- c(parts, sprintf("%d low-support predictions", n_low))
  }
  if (n_failed) {
    parts <- c(parts, sprintf("%d weight-normalization failures", n_failed))
  }
  if (n_incomplete) {
    parts <- c(parts, sprintf("%d incomplete prediction rows", n_incomplete))
  }
  if (length(parts)) {
    warning(
      paste0("EQuA prediction diagnostics: ", paste(parts, collapse = "; "), "."),
      call. = FALSE
    )
  }
}

#' Predict EQuA profiles
#'
#' @param object An `equa_fit` object.
#' @param new_estimated_age Numeric vector of estimated ages for new
#'   observations.
#' @param new_covariates Optional numeric covariates required by Kernel
#'   Ex-EQuA.
#' @param chronological_age Optional chronological ages. When supplied, qEAA
#'   profiles are returned in addition to qEA profiles.
#' @param support_action Optional override of the fitted support behavior.
#' @param min_effective_n Optional override of the low-support threshold.
#' @param nearest_k Optional override of the explicit nearest-neighbor fallback
#'   size.
#' @param ... Reserved for future methods.
#'
#' @return An object of class `equa_profile` containing qEA, optional qEAA,
#'   and support diagnostics.
#'
#' @examples
#' calibration <- subset(
#'   simulated_equa_data,
#'   sample_role == "calibration"
#' )
#' new_data <- subset(
#'   simulated_equa_data,
#'   sample_role == "prediction"
#' )
#' fit <- fit_equa(
#'   chronological_age = calibration$chronological_age,
#'   estimated_age = calibration$estimated_age,
#'   method = "kernel_end",
#'   probs = c(0.25, 0.5, 0.75)
#' )
#' profile <- predict(
#'   fit,
#'   new_estimated_age = new_data$estimated_age,
#'   chronological_age = new_data$chronological_age
#' )
#' head(as.data.frame(profile))
#'
#' @export
predict.equa_fit <- function(
    object,
    new_estimated_age,
    new_covariates = NULL,
    chronological_age = NULL,
    support_action = object$support_action,
    min_effective_n = object$min_effective_n,
    nearest_k = object$nearest_k,
    ...) {
  if (!inherits(object, "equa_fit")) {
    stop("`object` must inherit from `equa_fit`.", call. = FALSE)
  }
  support <- .validate_support_settings(
    support_action,
    min_effective_n,
    nearest_k
  )
  prediction <- .prepare_prediction_conditioning(
    object,
    new_estimated_age,
    new_covariates
  )
  n <- length(prediction$estimated_age)

  if (!is.null(chronological_age)) {
    chronological_age <- .validate_numeric_vector(
      chronological_age,
      "chronological_age",
      expected_length = n
    )
  }

  n_probs <- length(object$probs)
  qead <- matrix(NA_real_, nrow = n, ncol = n_probs)
  effective_sample_size <- rep(NA_real_, n)
  maximum_normalized_weight <- rep(NA_real_, n)
  nearest_standardized_distance <- rep(NA_real_, n)
  normalization_failed <- rep(FALSE, n)
  fallback_used <- rep(FALSE, n)

  outside_by_variable <- matrix(
    FALSE,
    nrow = n,
    ncol = length(object$conditioning_names),
    dimnames = list(NULL, object$conditioning_names)
  )
  if (ncol(outside_by_variable)) {
    for (j in seq_len(ncol(outside_by_variable))) {
      outside_by_variable[, j] <-
        prediction$raw[, j] < object$conditioning_range["minimum", j] |
        prediction$raw[, j] > object$conditioning_range["maximum", j]
      outside_by_variable[!is.finite(prediction$raw[, j]), j] <- NA
    }
  }

  valid_rows <- which(prediction$valid)
  if (object$method == "po") {
    weights <- rep(1 / object$n_calibration, object$n_calibration)
    unconditional <- .weighted_quantile_left(
      object$calibration_ead,
      weights,
      object$probs
    )
    qead[valid_rows, ] <- matrix(
      unconditional,
      nrow = length(valid_rows),
      ncol = n_probs,
      byrow = TRUE
    )
    effective_sample_size[valid_rows] <- object$n_calibration
    maximum_normalized_weight[valid_rows] <- 1 / object$n_calibration
  } else {
    for (i in valid_rows) {
      result <- .kernel_prediction_row(object, prediction$scaled[i, ])
      effective_sample_size[[i]] <- result$effective_sample_size
      maximum_normalized_weight[[i]] <- result$maximum_normalized_weight
      nearest_standardized_distance[[i]] <-
        result$nearest_standardized_distance
      normalization_failed[[i]] <- result$normalization_failed

      low <- is.finite(result$effective_sample_size) &&
        result$effective_sample_size < support$min_effective_n
      use_nearest <- support$support_action == "nearest" &&
        (isTRUE(low) || result$normalization_failed)

      if (use_nearest) {
        distances <- result$euclidean_distances
        if (!any(is.finite(distances))) {
          distances <- sqrt(rowSums(
            sweep(
              object$calibration_conditioning,
              2L,
              prediction$scaled[i, ],
              "-"
            )^2
          ))
        }
        nearest <- order(distances)[
          seq_len(min(support$nearest_k, object$n_calibration))
        ]
        qead[i, ] <- .weighted_quantile_left(
          object$calibration_ead[nearest],
          rep(1, length(nearest)),
          object$probs
        )
        fallback_used[[i]] <- TRUE
      } else {
        qead[i, ] <- result$quantiles
      }
    }
  }

  low_support <- is.finite(effective_sample_size) &
    effective_sample_size < support$min_effective_n
  outside_training_range <- if (ncol(outside_by_variable)) {
    apply(outside_by_variable, 1L, function(x) any(x, na.rm = TRUE))
  } else {
    rep(FALSE, n)
  }

  diagnostics <- data.frame(
    outside_training_range = outside_training_range,
    effective_sample_size = effective_sample_size,
    nearest_standardized_distance = nearest_standardized_distance,
    maximum_normalized_weight = maximum_normalized_weight,
    low_support = low_support,
    weight_normalization_failed = normalization_failed,
    fallback_used = fallback_used,
    incomplete_conditioning = !prediction$valid,
    stringsAsFactors = FALSE
  )

  problematic <- low_support | normalization_failed
  if (support$support_action == "na") {
    qead[problematic, ] <- NA_real_
  } else if (support$support_action == "error" && any(problematic)) {
    stop(
      sprintf(
        "EQuA prediction stopped: %d rows had low or invalid kernel support.",
        sum(problematic)
      ),
      call. = FALSE
    )
  } else if (support$support_action == "warn") {
    .summarize_support_warning(diagnostics)
  }

  qea <- sweep(qead, 1L, prediction$estimated_age, "+")
  colnames(qea) <- .profile_column_names("qea", object$probs)
  qeaa <- NULL
  if (!is.null(chronological_age)) {
    qeaa <- sweep(qea, 1L, chronological_age, "-")
    qeaa[!is.finite(chronological_age), ] <- NA_real_
    colnames(qeaa) <- .profile_column_names("qeaa", object$probs)
  }

  profile <- list(
    method = object$method,
    probs = object$probs,
    estimated_age = prediction$estimated_age,
    chronological_age = chronological_age,
    qea = qea,
    qeaa = qeaa,
    diagnostics = diagnostics,
    outside_range_by_variable = outside_by_variable,
    support_action = support$support_action,
    min_effective_n = support$min_effective_n,
    call = match.call()
  )
  class(profile) <- "equa_profile"
  profile
}
