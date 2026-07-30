test_that("automatic bandwidth matches the reproducibility rule", {
  x <- c(1, 2, 4, 8, 16, 32)

  expect_equal(
    EQuA:::.select_bandwidth_rule(x, 1L),
    reference_bandwidth(x, 1L)
  )
  expect_equal(
    EQuA:::.select_bandwidth_rule(x, 3L),
    reference_bandwidth(x, 3L)
  )
})

test_that("bandwidth selection is fold-local and expressed in raw units", {
  age <- seq(20, 80, length.out = 30)
  estimated <- age + sin(age)
  z <- seq(-3, 6, length.out = 30)^3
  fit <- fit_equa(
    age,
    estimated,
    covariates = data.frame(exposure = z),
    method = "kernel_ex",
    min_effective_n = 1
  )

  expect_equal(
    unname(fit$bandwidth),
    c(
      reference_bandwidth(estimated, 2L),
      reference_bandwidth(z, 2L)
    )
  )
  expect_equal(
    unname(fit$bandwidth_scaled),
    unname(fit$bandwidth / fit$conditioning_scale)
  )
  expect_true(fit$bandwidth_selected)
})

test_that("kernel estimates are invariant to a change of covariate units", {
  age <- seq(18, 78, length.out = 80)
  estimated <- age + 2 * sin(age / 9)
  z <- seq(-2, 3, length.out = 80) + cos(age / 7)

  fit_original <- fit_equa(
    age,
    estimated,
    covariates = data.frame(z = z),
    method = "kernel_ex",
    probs = c(0.2, 0.5, 0.8),
    min_effective_n = 1
  )
  fit_rescaled <- fit_equa(
    age,
    estimated,
    covariates = data.frame(z = -1000 * z + 17),
    method = "kernel_ex",
    probs = c(0.2, 0.5, 0.8),
    min_effective_n = 1
  )

  prediction_original <- predict(
    fit_original,
    new_estimated_age = estimated[c(7, 29, 65)],
    new_covariates = data.frame(z = z[c(7, 29, 65)])
  )
  prediction_rescaled <- predict(
    fit_rescaled,
    new_estimated_age = estimated[c(7, 29, 65)],
    new_covariates = data.frame(z = -1000 * z[c(7, 29, 65)] + 17)
  )

  expect_equal(prediction_rescaled$qea, prediction_original$qea)
  expect_equal(
    unname(fit_rescaled$bandwidth[["z"]]),
    1000 * unname(fit_original$bandwidth[["z"]])
  )
})

test_that("age-unit changes transform profiles equivariantly", {
  age <- seq(18, 78, length.out = 80)
  estimated <- age + 2 * sin(age / 9)
  index <- c(7, 29, 65)

  fit_years <- fit_equa(
    age,
    estimated,
    method = "kernel_end",
    probs = c(0.2, 0.5, 0.8),
    min_effective_n = 1
  )
  fit_months <- fit_equa(
    12 * age,
    12 * estimated,
    method = "kernel_end",
    probs = c(0.2, 0.5, 0.8),
    min_effective_n = 1
  )

  profile_years <- predict(
    fit_years,
    new_estimated_age = estimated[index],
    chronological_age = age[index]
  )
  profile_months <- predict(
    fit_months,
    new_estimated_age = 12 * estimated[index],
    chronological_age = 12 * age[index]
  )

  expect_equal(profile_months$qea, 12 * profile_years$qea)
  expect_equal(profile_months$qeaa, 12 * profile_years$qeaa)
})
