test_that("Po-EQuA matches a hand-computed empirical quantile", {
  fit <- fit_equa(
    chronological_age = c(10, 13, 16, 19),
    estimated_age = c(9, 10, 13, 17),
    method = "po",
    probs = c(0.25, 0.5, 0.75),
    min_effective_n = 1
  )
  profile <- predict(
    fit,
    new_estimated_age = c(20, 30),
    chronological_age = c(22, 31)
  )

  expect_equal(unname(profile$qea), rbind(c(21, 22, 23), c(31, 32, 33)))
  expect_equal(unname(profile$qeaa), rbind(c(-1, 0, 1), c(0, 1, 2)))
  expect_identical(
    colnames(profile$qea),
    c("qea_p25", "qea_p50", "qea_p75")
  )
  expect_equal(profile$diagnostics$effective_sample_size, c(4, 4))
  expect_false(any(profile$diagnostics$fallback_used))
})

test_that("Kernel End-EQuA matches the frozen Gaussian-kernel calculation", {
  age <- c(20, 25, 31, 38, 46, 55)
  estimated <- c(19, 27, 30, 40, 43, 58)
  ead <- age - estimated
  probs <- c(0.2, 0.5, 0.8)
  bandwidth <- 8
  new_estimated <- c(29, 47)

  fit <- fit_equa(
    age,
    estimated,
    method = "kernel_end",
    probs = probs,
    bandwidth = bandwidth,
    min_effective_n = 1
  )
  profile <- predict(fit, new_estimated_age = new_estimated)

  expected_qead <- t(vapply(
    new_estimated,
    function(value) {
      weights <- exp(-0.5 * ((estimated - value) / bandwidth)^2)
      reference_weighted_quantile(ead, weights, probs)
    },
    numeric(length(probs))
  ))
  expect_equal(unname(profile$qea), expected_qead + new_estimated)
  expect_false(any(profile$diagnostics$fallback_used))
})

test_that("Kernel Ex-EQuA uses the product Gaussian kernel", {
  age <- c(20, 25, 31, 38, 46, 55)
  estimated <- c(19, 27, 30, 40, 43, 58)
  z <- c(-1.2, 0.1, -0.4, 1.1, 0.7, 1.8)
  ead <- age - estimated
  probs <- c(0.25, 0.75)
  bandwidth <- c(9, 0.8)
  new_estimated <- 35
  new_z <- 0.25

  fit <- fit_equa(
    age,
    estimated,
    covariates = data.frame(z = z),
    method = "kernel_ex",
    probs = probs,
    bandwidth = bandwidth,
    min_effective_n = 1
  )
  profile <- predict(
    fit,
    new_estimated_age = new_estimated,
    new_covariates = data.frame(z = new_z)
  )

  weights <- exp(-0.5 * (
    ((estimated - new_estimated) / bandwidth[[1L]])^2 +
      ((z - new_z) / bandwidth[[2L]])^2
  ))
  expected <- new_estimated +
    reference_weighted_quantile(ead, weights, probs)
  expect_equal(as.numeric(profile$qea), expected)
  expect_false(is.unsorted(as.numeric(profile$qea)))
  expect_named(
    profile$diagnostics,
    c(
      "outside_training_range",
      "effective_sample_size",
      "nearest_standardized_distance",
      "maximum_normalized_weight",
      "low_support",
      "weight_normalization_failed",
      "fallback_used",
      "incomplete_conditioning"
    )
  )
})

test_that("low support never triggers an implicit nearest-neighbor fallback", {
  age <- seq(20, 60, length.out = 20)
  estimated <- age + sin(age)
  fit <- fit_equa(
    age,
    estimated,
    method = "kernel_end",
    bandwidth = 0.01,
    min_effective_n = 10,
    support_action = "warn"
  )

  profile <- NULL
  expect_warning(
    profile <- predict(fit, new_estimated_age = 1000),
    "low-support"
  )
  expect_true(profile$diagnostics$low_support)
  expect_false(profile$diagnostics$fallback_used)
  expect_true(all(is.finite(profile$qea)))

  profile_na <- predict(
    fit,
    new_estimated_age = 1000,
    support_action = "na"
  )
  expect_true(all(is.na(profile_na$qea)))
  expect_error(
    predict(fit, new_estimated_age = 1000, support_action = "error"),
    "low or invalid kernel support"
  )

  profile_nearest <- predict(
    fit,
    new_estimated_age = 1000,
    support_action = "nearest",
    nearest_k = 5
  )
  expect_true(profile_nearest$diagnostics$fallback_used)
  expect_true(all(is.finite(profile_nearest$qea)))
})

test_that("range, missing-value, and covariate diagnostics are explicit", {
  age <- seq(20, 70, length.out = 30)
  estimated <- age + cos(age)
  z <- seq(-1, 1, length.out = 30)
  fit <- fit_equa(
    age,
    estimated,
    covariates = data.frame(z = z),
    method = "kernel_ex",
    min_effective_n = 1
  )

  profile <- NULL
  expect_warning(
    profile <- predict(
      fit,
      new_estimated_age = c(estimated[[1L]] - 10, NA_real_),
      new_covariates = data.frame(z = c(4, 0))
    ),
    "incomplete prediction rows"
  )
  expect_true(profile$diagnostics$outside_training_range[[1L]])
  expect_true(all(profile$outside_range_by_variable[1, ]))
  expect_true(profile$diagnostics$incomplete_conditioning[[2L]])
  expect_true(all(is.na(profile$qea[2, ])))
  expect_error(
    predict(
      fit,
      new_estimated_age = 30,
      new_covariates = data.frame(wrong_name = 0)
    ),
    "column names and order"
  )
})

test_that("calibration omission is method-specific and documented in the fit", {
  age <- c(20, 30, NA, 50, 60)
  estimated <- c(21, 29, 40, Inf, 61)

  fit <- NULL
  expect_warning(
    fit <- fit_equa(
      age,
      estimated,
      method = "po",
      min_effective_n = 1
    ),
    "Omitted 2"
  )
  expect_equal(fit$n_input, 5)
  expect_equal(fit$n_calibration, 3)
  expect_equal(fit$n_omitted, 2)
  expect_error(
    fit_equa(age, estimated, method = "po", na_action = "fail"),
    "2 incomplete"
  )
})

test_that("non-finite chronological age does not contaminate qEA", {
  fit <- fit_equa(
    chronological_age = c(10, 13, 16, 19),
    estimated_age = c(9, 10, 13, 17),
    method = "po",
    probs = c(0.25, 0.75),
    min_effective_n = 1
  )
  profile <- predict(
    fit,
    new_estimated_age = c(20, 30),
    chronological_age = c(NA_real_, Inf)
  )

  expect_true(all(is.finite(profile$qea)))
  expect_true(all(is.na(profile$qeaa)))
})
