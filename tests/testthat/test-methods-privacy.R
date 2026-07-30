test_that("fit print and summary methods do not reveal participant values", {
  age <- c(10001, 10003, 10008, 10012, 10020)
  estimated <- c(9999, 10002, 10005, 10015, 10018)
  fit <- fit_equa(
    age,
    estimated,
    method = "po",
    probs = c(0.25, 0.75),
    min_effective_n = 1
  )

  printed <- capture.output(print(fit))
  summarized <- capture.output(print(summary(fit)))
  expect_match(paste(printed, collapse = " "), "Privacy")
  expect_match(paste(summarized, collapse = " "), "Privacy")
  expect_false(grepl("10001|10003|10008|10012|10020", paste(printed, collapse = " ")))
  expect_false(grepl(
    "10001|10003|10008|10012|10020",
    paste(summarized, collapse = " ")
  ))
  expect_false("subject_id" %in% names(fit))
})

test_that("profile coercion supports wide and long layouts", {
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
    chronological_age = c(22, 31)
  )

  wide <- as.data.frame(profile)
  long <- as.data.frame(profile, format = "long")
  expect_equal(nrow(wide), 2)
  expect_equal(nrow(long), 4)
  expect_named(wide, c(
    "observation",
    "estimated_age",
    "chronological_age",
    "qea_p25",
    "qea_p75",
    "qeaa_p25",
    "qeaa_p75",
    names(profile$diagnostics)
  ))
  expect_setequal(unique(long$tau), c(0.25, 0.75))
})

test_that("profile plots are generated with base graphics", {
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
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot(profile, type = "qea"))
  expect_true(file.exists(output))
})

test_that("cross-fit objects reuse profile methods through inheritance", {
  age <- seq(20, 70, length.out = 24)
  estimated <- age + sin(age / 4)
  profile <- crossfit_equa(
    age,
    estimated,
    method = "po",
    probs = c(0.25, 0.5, 0.75),
    folds = 4,
    repeats = 1,
    seed = 52,
    min_effective_n = 1
  )

  printed <- capture.output(print(profile))
  summarized <- summary(profile)
  long <- as.data.frame(profile, format = "long")
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_match(paste(printed, collapse = " "), "Folds:")
  expect_identical(summarized$config$folds, 4L)
  expect_equal(nrow(long), 24 * 3)
  expect_invisible(plot(profile, type = "qeaa", rows = 1:2))
})

test_that("diagnostics can include per-variable range flags", {
  age <- seq(20, 60, length.out = 20)
  estimated <- age + sin(age)
  fit <- fit_equa(
    age,
    estimated,
    method = "kernel_end",
    min_effective_n = 1
  )
  profile <- predict(fit, new_estimated_age = c(30, 90))
  diagnostics <- equa_diagnostics(profile, by_variable = TRUE)

  expect_named(
    diagnostics,
    c("summary", "outside_range_by_variable")
  )
  expect_identical(
    colnames(diagnostics$outside_range_by_variable),
    "estimated_age"
  )
  expect_true(diagnostics$outside_range_by_variable[2, 1])
})
