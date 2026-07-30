reference_end_crossfit <- function(
    age,
    estimated,
    probs,
    folds,
    repeats,
    seed) {
  ead <- age - estimated
  runs <- array(
    NA_real_,
    dim = c(length(age), length(probs), repeats)
  )
  assignments <- matrix(NA_integer_, nrow = length(age), ncol = repeats)

  for (repeat_id in seq_len(repeats)) {
    assignment <- reference_age_folds(
      age,
      folds,
      seed + 100L * repeat_id
    )
    assignments[, repeat_id] <- assignment
    for (fold_id in seq_len(folds)) {
      test <- which(assignment == fold_id)
      train <- which(assignment != fold_id)
      bandwidth <- reference_bandwidth(estimated[train], 1L)
      for (i in test) {
        weights <- exp(
          -0.5 * ((estimated[train] - estimated[[i]]) / bandwidth)^2
        )
        runs[i, , repeat_id] <- estimated[[i]] +
          reference_weighted_quantile(ead[train], weights, probs)
      }
    }
  }

  list(
    qea = apply(runs, c(1L, 2L), stats::median),
    assignments = assignments
  )
}

test_that("fold construction matches the frozen stratified algorithm", {
  age <- c(50, 20, 30, 30, 80, 42, 60, 26, 70, 38, 48, 55)

  set.seed(313)
  before <- .Random.seed
  result <- EQuA:::.make_age_stratified_folds(age, 3, 909, 10)
  expect_identical(.Random.seed, before)
  expect_identical(result, reference_age_folds(age, 3, 909, 10))
  expect_setequal(unique(result), 1:3)
})

test_that("repeated End-EQuA cross-fitting matches a reference implementation", {
  age <- seq(18, 82, length.out = 60)
  estimated <- age + 2.5 * sin(age / 8) + seq(-0.5, 0.5, length.out = 60)
  probs <- c(0.2, 0.5, 0.8)
  reference <- reference_end_crossfit(
    age,
    estimated,
    probs,
    folds = 5,
    repeats = 2,
    seed = 2026
  )

  set.seed(77)
  state_before <- .Random.seed
  result <- crossfit_equa(
    age,
    estimated,
    method = "kernel_end",
    probs = probs,
    folds = 5,
    repeats = 2,
    seed = 2026,
    min_effective_n = 1,
    keep_repeats = TRUE
  )

  expect_identical(.Random.seed, state_before)
  expect_s3_class(result, "equa_crossfit")
  expect_s3_class(result, "equa_profile")
  expect_equal(unname(result$qea), reference$qea)
  expect_identical(unname(result$fold_assignments), reference$assignments)
  expect_equal(
    unname(result$qeaa),
    unname(sweep(result$qea, 1L, age, "-"))
  )
  expect_equal(dim(result$repeated_profiles$qea), c(60, 3, 2))
  expect_equal(nrow(result$fold_diagnostics), 10)
  expect_equal(nrow(result$bandwidth_diagnostics), 10)
  expect_identical(result$config$rng_kind, RNGkind())
  expect_identical(result$config$r_version, as.character(getRversion()))
})

test_that("cross-fitting can omit large audit components", {
  age <- seq(20, 70, length.out = 40)
  estimated <- age + sin(age / 5)

  result <- crossfit_equa(
    age,
    estimated,
    method = "po",
    probs = c(0.25, 0.75),
    folds = 4,
    repeats = 2,
    seed = 11,
    min_effective_n = 1,
    keep_repeats = FALSE,
    keep_fold_assignments = FALSE
  )

  expect_null(result$repeated_profiles)
  expect_null(result$fold_assignments)
  expect_equal(nrow(result$bandwidth_diagnostics), 0)
  expect_equal(nrow(result$diagnostics), length(age))
})

test_that("Ex-EQuA cross-fitting handles multiple conditioning variables", {
  age <- seq(20, 80, length.out = 60)
  estimated <- age + 2 * sin(age / 8)
  covariates <- data.frame(
    z1 = seq(-2, 2, length.out = 60),
    z2 = cos(age / 6)
  )

  result <- crossfit_equa(
    age,
    estimated,
    covariates = covariates,
    method = "kernel_ex",
    probs = c(0.2, 0.5, 0.8),
    folds = 5,
    repeats = 1,
    seed = 87,
    min_effective_n = 1
  )

  expect_equal(dim(result$qea), c(60, 3))
  expect_true(all(apply(result$qea, 1L, function(x) !is.unsorted(x))))
  expect_identical(
    colnames(result$outside_range_by_variable),
    c("estimated_age", "z1", "z2")
  )
  expect_equal(nrow(result$bandwidth_diagnostics), 15)
})

test_that("cross-fitting restores omitted rows as explicit missing profiles", {
  age <- c(seq(20, 70, length.out = 40), NA_real_)
  estimated <- c(age[1:40] + cos(age[1:40]), 45)

  result <- NULL
  expect_warning(
    result <- crossfit_equa(
      age,
      estimated,
      method = "kernel_end",
      probs = c(0.25, 0.75),
      folds = 4,
      repeats = 2,
      seed = 22,
      min_effective_n = 1
    ),
    "Omitted 1"
  )

  expect_equal(nrow(result$qea), 41)
  expect_true(all(is.na(result$qea[41, ])))
  expect_true(result$diagnostics$incomplete_conditioning[[41L]])
  expect_true(all(is.na(result$fold_assignments[41, ])))
})
