.with_preserved_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(code)
}

.make_age_stratified_folds <- function(age, folds, seed, n_strata) {
  n <- length(age)
  if (n < folds) {
    stop("The number of complete rows must be at least `folds`.", call. = FALSE)
  }

  # Match the paper workflow for ordinary cohort sizes: up to ten
  # approximately equal chronological-age strata, randomized within stratum,
  # then cyclic fold assignment. For small inputs, reduce the number of strata
  # so every stratum can populate all folds.
  n_strata_use <- min(n_strata, max(1L, floor(n / folds)))
  age_order <- order(age, seq_along(age))
  stratum <- integer(n)
  # This is the integer form of dplyr::ntile() used by the paper workflow:
  # larger strata occur first when n is not divisible by n_strata_use.
  stratum[age_order] <-
    ((seq_len(n) - 1) * as.double(n_strata_use)) %/% n + 1L

  .with_preserved_seed(seed, {
    assignment <- integer(n)
    for (s in seq_len(n_strata_use)) {
      # Preserve chronological-age order inside each stratum before applying
      # the random permutation, matching the arranged input in the original
      # dplyr workflow.
      index <- age_order[stratum[age_order] == s]
      shuffled <- sample(index, length(index))
      assignment[shuffled] <- rep(seq_len(folds), length.out = length(index))
    }
    assignment
  })
}

.aggregate_profile_array <- function(x, aggregate) {
  fun <- if (aggregate == "median") .median_or_na else .mean_or_na
  result <- apply(x, c(1L, 2L), fun)
  if (is.null(dim(result))) {
    result <- matrix(result, nrow = dim(x)[[1L]], ncol = dim(x)[[2L]])
  }
  result
}

#' Repeated cross-fitted EQuA profiles
#'
#' Produces strictly out-of-fold qEA and qEAA profiles using repeated,
#' chronological-age-stratified K-fold cross-fitting.
#'
#' @inheritParams fit_equa
#' @param folds Number of folds.
#' @param repeats Number of repeated fold assignments.
#' @param seed Integer base seed. Repeat `r` uses `seed + 100 * r`, matching the
#'   reproducibility workflow.
#' @param n_strata Maximum number of chronological-age strata.
#' @param aggregate Pointwise aggregation across repeats, either `"median"` or
#'   `"mean"`.
#' @param keep_repeats Whether to retain all repeat-level qEA and qEAA arrays.
#' @param keep_fold_assignments Whether to retain the fold-assignment matrix.
#'
#' @return An object with classes `equa_crossfit` and `equa_profile`.
#'   Cross-fitting results contain participant-level inputs and profiles; review
#'   disclosure risk before saving or sharing them.
#'
#' @examples
#' calibration <- subset(
#'   simulated_equa_data,
#'   sample_role == "calibration"
#' )
#' profile <- crossfit_equa(
#'   chronological_age = calibration$chronological_age,
#'   estimated_age = calibration$estimated_age,
#'   method = "kernel_end",
#'   probs = c(0.25, 0.5, 0.75),
#'   folds = 5,
#'   repeats = 2,
#'   seed = 2026
#' )
#' profile
#'
#' @export
crossfit_equa <- function(
    chronological_age,
    estimated_age,
    covariates = NULL,
    method = c("kernel_end", "po", "kernel_ex"),
    probs = seq(0.1, 0.9, by = 0.1),
    bandwidth = NULL,
    folds = 5,
    repeats = 50,
    seed = 2026,
    n_strata = 10,
    aggregate = c("median", "mean"),
    na_action = c("omit", "fail"),
    support_action = c("warn", "na", "error", "nearest"),
    min_effective_n = 10,
    nearest_k = 30,
    keep_repeats = FALSE,
    keep_fold_assignments = TRUE) {
  method <- .match_equa_method(method)
  probs <- .validate_probs(probs)
  aggregate <- match.arg(aggregate)
  na_action <- match.arg(na_action)
  support <- .validate_support_settings(
    support_action,
    min_effective_n,
    nearest_k
  )

  integer_settings <- list(
    folds = folds,
    repeats = repeats,
    seed = seed,
    n_strata = n_strata
  )
  for (name in names(integer_settings)) {
    value <- integer_settings[[name]]
    if (!is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value) ||
        value != floor(value) ||
        value > .Machine$integer.max) {
      stop(sprintf("`%s` must be one finite integer.", name), call. = FALSE)
    }
  }
  folds <- as.integer(folds)
  repeats <- as.integer(repeats)
  seed <- as.integer(seed)
  n_strata <- as.integer(n_strata)
  if (folds < 2L) stop("`folds` must be at least 2.", call. = FALSE)
  if (repeats < 1L) stop("`repeats` must be at least 1.", call. = FALSE)
  if (seed < 0L) stop("`seed` must be non-negative.", call. = FALSE)
  if (n_strata < 1L) stop("`n_strata` must be at least 1.", call. = FALSE)
  if (as.double(seed) + 100 * as.double(repeats) >
      .Machine$integer.max) {
    stop("`seed` is too large for the requested number of repeats.", call. = FALSE)
  }
  if (!is.logical(keep_repeats) || length(keep_repeats) != 1L || is.na(keep_repeats)) {
    stop("`keep_repeats` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(keep_fold_assignments) ||
      length(keep_fold_assignments) != 1L ||
      is.na(keep_fold_assignments)) {
    stop("`keep_fold_assignments` must be TRUE or FALSE.", call. = FALSE)
  }

  inputs <- .prepare_calibration_inputs(
    chronological_age,
    estimated_age,
    covariates,
    method,
    na_action
  )
  n_total <- inputs$n_input
  valid_index <- which(inputs$complete)
  n_valid <- length(valid_index)
  if (n_valid < folds + 2L) {
    stop(
      "Complete rows must exceed the number of folds by at least two.",
      call. = FALSE
    )
  }

  n_probs <- length(probs)
  conditioning_dimension <- if (method == "po") {
    0L
  } else if (method == "kernel_end") {
    1L
  } else {
    1L + ncol(inputs$covariates)
  }

  qea_runs <- array(NA_real_, dim = c(n_valid, n_probs, repeats))
  qeaa_runs <- array(NA_real_, dim = c(n_valid, n_probs, repeats))
  ess_runs <- matrix(NA_real_, nrow = n_valid, ncol = repeats)
  distance_runs <- matrix(NA_real_, nrow = n_valid, ncol = repeats)
  max_weight_runs <- matrix(NA_real_, nrow = n_valid, ncol = repeats)
  low_runs <- matrix(FALSE, nrow = n_valid, ncol = repeats)
  failure_runs <- matrix(FALSE, nrow = n_valid, ncol = repeats)
  fallback_runs <- matrix(FALSE, nrow = n_valid, ncol = repeats)
  outside_runs <- if (conditioning_dimension) {
    array(FALSE, dim = c(n_valid, conditioning_dimension, repeats))
  } else {
    NULL
  }
  fold_assignments_valid <- matrix(
    NA_integer_,
    nrow = n_valid,
    ncol = repeats
  )

  bandwidth_rows <- list()
  fold_rows <- list()
  bandwidth_at <- 0L
  fold_at <- 0L
  rng_kind <- RNGkind()
  r_version <- as.character(getRversion())

  for (repeat_id in seq_len(repeats)) {
    repeat_seed <- as.integer(as.double(seed) + 100 * repeat_id)
    assignment <- .make_age_stratified_folds(
      inputs$chronological_age,
      folds,
      repeat_seed,
      n_strata
    )
    fold_assignments_valid[, repeat_id] <- assignment

    for (fold_id in seq_len(folds)) {
      test_index <- which(assignment == fold_id)
      train_index <- which(assignment != fold_id)
      train_covariates <- if (is.null(inputs$covariates)) {
        NULL
      } else {
        inputs$covariates[train_index, , drop = FALSE]
      }
      test_covariates <- if (is.null(inputs$covariates)) {
        NULL
      } else {
        inputs$covariates[test_index, , drop = FALSE]
      }

      fitted <- fit_equa(
        chronological_age = inputs$chronological_age[train_index],
        estimated_age = inputs$estimated_age[train_index],
        covariates = train_covariates,
        method = method,
        probs = probs,
        bandwidth = bandwidth,
        na_action = "fail",
        support_action = support$support_action,
        min_effective_n = support$min_effective_n,
        nearest_k = support$nearest_k
      )

      predicted <- if (support$support_action == "warn") {
        suppressWarnings(predict(
          fitted,
          new_estimated_age = inputs$estimated_age[test_index],
          new_covariates = test_covariates,
          chronological_age = inputs$chronological_age[test_index]
        ))
      } else {
        predict(
          fitted,
          new_estimated_age = inputs$estimated_age[test_index],
          new_covariates = test_covariates,
          chronological_age = inputs$chronological_age[test_index]
        )
      }

      qea_runs[test_index, , repeat_id] <- predicted$qea
      qeaa_runs[test_index, , repeat_id] <- predicted$qeaa
      ess_runs[test_index, repeat_id] <-
        predicted$diagnostics$effective_sample_size
      distance_runs[test_index, repeat_id] <-
        predicted$diagnostics$nearest_standardized_distance
      max_weight_runs[test_index, repeat_id] <-
        predicted$diagnostics$maximum_normalized_weight
      low_runs[test_index, repeat_id] <- predicted$diagnostics$low_support
      failure_runs[test_index, repeat_id] <-
        predicted$diagnostics$weight_normalization_failed
      fallback_runs[test_index, repeat_id] <-
        predicted$diagnostics$fallback_used
      if (conditioning_dimension) {
        outside_runs[test_index, , repeat_id] <-
          predicted$outside_range_by_variable
      }

      if (length(fitted$bandwidth)) {
        for (j in seq_along(fitted$bandwidth)) {
          bandwidth_at <- bandwidth_at + 1L
          bandwidth_rows[[bandwidth_at]] <- data.frame(
            repeat_id = repeat_id,
            fold = fold_id,
            conditioning_variable = names(fitted$bandwidth)[[j]],
            bandwidth = unname(fitted$bandwidth[[j]]),
            training_n = length(train_index),
            test_n = length(test_index),
            stringsAsFactors = FALSE
          )
        }
      }

      fold_at <- fold_at + 1L
      fold_rows[[fold_at]] <- data.frame(
        repeat_id = repeat_id,
        fold = fold_id,
        training_n = length(train_index),
        test_n = length(test_index),
        n_low_support = sum(predicted$diagnostics$low_support),
        n_weight_failures = sum(
          predicted$diagnostics$weight_normalization_failed
        ),
        n_fallback = sum(predicted$diagnostics$fallback_used),
        median_effective_sample_size = .median_or_na(
          predicted$diagnostics$effective_sample_size
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  if (support$support_action == "warn") {
    n_low <- sum(low_runs)
    n_failed <- sum(failure_runs)
    if (n_low || n_failed) {
      warning(
        sprintf(
          paste(
            "Cross-fitting diagnostics across all repeats:",
            "%d low-support predictions and %d weight-normalization failures."
          ),
          n_low,
          n_failed
        ),
        call. = FALSE
      )
    }
  }

  qea_valid <- .aggregate_profile_array(qea_runs, aggregate)
  qeaa_valid <- .aggregate_profile_array(qeaa_runs, aggregate)
  qea <- matrix(NA_real_, nrow = n_total, ncol = n_probs)
  qeaa <- matrix(NA_real_, nrow = n_total, ncol = n_probs)
  qea[valid_index, ] <- qea_valid
  qeaa[valid_index, ] <- qeaa_valid
  colnames(qea) <- .profile_column_names("qea", probs)
  colnames(qeaa) <- .profile_column_names("qeaa", probs)

  diagnostics <- data.frame(
    outside_training_range = rep(FALSE, n_total),
    effective_sample_size = rep(NA_real_, n_total),
    nearest_standardized_distance = rep(NA_real_, n_total),
    maximum_normalized_weight = rep(NA_real_, n_total),
    low_support = rep(FALSE, n_total),
    low_support_rate = rep(NA_real_, n_total),
    weight_normalization_failed = rep(FALSE, n_total),
    fallback_used = rep(FALSE, n_total),
    incomplete_conditioning = !inputs$complete,
    stringsAsFactors = FALSE
  )
  diagnostics$effective_sample_size[valid_index] <-
    apply(ess_runs, 1L, .median_or_na)
  diagnostics$nearest_standardized_distance[valid_index] <-
    apply(distance_runs, 1L, .median_or_na)
  diagnostics$maximum_normalized_weight[valid_index] <-
    apply(max_weight_runs, 1L, .median_or_na)
  diagnostics$low_support[valid_index] <- apply(low_runs, 1L, any)
  diagnostics$low_support_rate[valid_index] <- rowMeans(low_runs)
  diagnostics$weight_normalization_failed[valid_index] <-
    apply(failure_runs, 1L, any)
  diagnostics$fallback_used[valid_index] <- apply(fallback_runs, 1L, any)

  conditioning_names <- if (conditioning_dimension) {
    if (method == "kernel_end") {
      "estimated_age"
    } else {
      c("estimated_age", colnames(inputs$covariates))
    }
  } else {
    character()
  }
  outside_by_variable <- matrix(
    FALSE,
    nrow = n_total,
    ncol = conditioning_dimension,
    dimnames = list(NULL, conditioning_names)
  )
  if (conditioning_dimension) {
    outside_valid <- apply(outside_runs, c(1L, 2L), any)
    if (is.null(dim(outside_valid))) {
      outside_valid <- matrix(
        outside_valid,
        nrow = n_valid,
        ncol = conditioning_dimension
      )
    }
    outside_by_variable[valid_index, ] <- outside_valid
    diagnostics$outside_training_range[valid_index] <-
      apply(outside_valid, 1L, any)
  }

  full_fold_assignments <- NULL
  if (keep_fold_assignments) {
    full_fold_assignments <- matrix(
      NA_integer_,
      nrow = n_total,
      ncol = repeats,
      dimnames = list(NULL, paste0("repeat_", seq_len(repeats)))
    )
    full_fold_assignments[valid_index, ] <- fold_assignments_valid
  }

  repeated_profiles <- NULL
  if (keep_repeats) {
    qea_full <- array(
      NA_real_,
      dim = c(n_total, n_probs, repeats),
      dimnames = list(
        NULL,
        .profile_column_names("qea", probs),
        paste0("repeat_", seq_len(repeats))
      )
    )
    qeaa_full <- array(
      NA_real_,
      dim = c(n_total, n_probs, repeats),
      dimnames = list(
        NULL,
        .profile_column_names("qeaa", probs),
        paste0("repeat_", seq_len(repeats))
      )
    )
    qea_full[valid_index, , ] <- qea_runs
    qeaa_full[valid_index, , ] <- qeaa_runs
    repeated_profiles <- list(qea = qea_full, qeaa = qeaa_full)
  }

  bandwidth_diagnostics <- if (length(bandwidth_rows)) {
    do.call(rbind, bandwidth_rows)
  } else {
    data.frame(
      repeat_id = integer(),
      fold = integer(),
      conditioning_variable = character(),
      bandwidth = numeric(),
      training_n = integer(),
      test_n = integer()
    )
  }

  profile <- list(
    method = method,
    probs = probs,
    estimated_age = .validate_numeric_vector(
      estimated_age,
      "estimated_age",
      expected_length = n_total
    ),
    chronological_age = .validate_numeric_vector(
      chronological_age,
      "chronological_age",
      expected_length = n_total
    ),
    qea = qea,
    qeaa = qeaa,
    diagnostics = diagnostics,
    outside_range_by_variable = outside_by_variable,
    support_action = support$support_action,
    min_effective_n = support$min_effective_n,
    n_input = n_total,
    n_complete = n_valid,
    n_omitted = inputs$n_omitted,
    config = list(
      folds = folds,
      repeats = repeats,
      seed = seed,
      n_strata = n_strata,
      aggregate = aggregate,
      bandwidth = bandwidth,
      na_action = na_action,
      support_action = support$support_action,
      min_effective_n = support$min_effective_n,
      nearest_k = support$nearest_k,
      rng_kind = rng_kind,
      r_version = r_version,
      keep_repeats = keep_repeats,
      keep_fold_assignments = keep_fold_assignments
    ),
    fold_assignments = full_fold_assignments,
    fold_diagnostics = do.call(rbind, fold_rows),
    bandwidth_diagnostics = bandwidth_diagnostics,
    repeated_profiles = repeated_profiles,
    call = match.call(),
    privacy_notice = paste(
      "This object contains participant-level input and profile values;",
      "review disclosure risk before saving or sharing it."
    )
  )
  class(profile) <- c("equa_crossfit", "equa_profile")
  profile
}
