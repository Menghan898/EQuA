test_that("fit inputs are validated before calculation", {
  expect_error(
    fit_equa(1:4, 1:3, method = "po"),
    "length 4"
  )
  expect_error(
    fit_equa(1:4, 1:4, method = "kernel_ex"),
    "covariates.*required"
  )
  expect_error(
    fit_equa(1:4, 1:4, covariates = 1:4, method = "po"),
    "must be NULL"
  )
  expect_error(
    fit_equa(1:4, 1:4, method = "po", probs = c(0.5, 0.25)),
    "strictly increasing"
  )
  expect_error(
    fit_equa(1:4, 1:4, method = "kernel_end", bandwidth = 0),
    "positive and finite"
  )
})

test_that("constant conditioning variables are rejected clearly", {
  expect_error(
    fit_equa(1:5, rep(3, 5), method = "kernel_end"),
    "positive variation.*estimated_age"
  )
  expect_error(
    fit_equa(
      1:5,
      1:5,
      covariates = data.frame(z = rep(1, 5)),
      method = "kernel_ex"
    ),
    "positive variation.*z"
  )
})

test_that("cross-fitting settings are validated", {
  age <- seq(20, 60, length.out = 20)
  estimated <- age + sin(age)

  expect_error(
    crossfit_equa(age, estimated, folds = 1),
    "at least 2"
  )
  expect_error(
    crossfit_equa(age, estimated, repeats = 0),
    "at least 1"
  )
  expect_error(
    crossfit_equa(age, estimated, folds = 19),
    "exceed the number of folds"
  )
  expect_error(
    crossfit_equa(
      age,
      estimated,
      repeats = 1,
      seed = .Machine$integer.max
    ),
    "seed.*too large"
  )
  expect_error(
    fit_equa(
      age,
      estimated,
      method = "kernel_end",
      nearest_k = .Machine$integer.max + 1
    ),
    "positive integer"
  )
})

test_that("profile column names stay unique for close probabilities", {
  probs <- c(0.5, 0.500000000000001)
  names <- EQuA:::.profile_column_names("qea", probs)

  expect_length(unique(names), 2)
})
