test_that("package metadata points to the fixed repository", {
  metadata <- utils::packageDescription("EQuA")

  expect_match(metadata$URL, "https://github.com/Menghan898/EQuA", fixed = TRUE)
  expect_match(
    metadata$BugReports,
    "https://github.com/Menghan898/EQuA/issues",
    fixed = TRUE
  )
  expect_match(
    metadata[["Authors@R"]],
    'email = "menghany@umich.edu"',
    fixed = TRUE
  )
  authors <- eval(parse(text = metadata[["Authors@R"]]))
  expect_length(authors, 1L)
  expect_identical(authors[[1L]]$given, "Menghan")
  expect_identical(authors[[1L]]$family, "Yi")
  expect_setequal(authors[[1L]]$role, c("aut", "cre"))
  expect_identical(metadata$License, "MIT + file LICENSE")
})

test_that("simulated example data are available and non-sensitive", {
  dat <- EQuA::simulated_equa_data

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 240L)
  expect_named(
    dat,
    c(
      "sample_id",
      "sample_role",
      "chronological_age",
      "estimated_age",
      "z"
    )
  )
  expect_setequal(unique(dat$sample_role), c("calibration", "prediction"))
})

test_that("public GEO examples have the documented compact interface", {
  examples <- list(
    gse40279_equa = EQuA::gse40279_equa,
    gse50660_equa = EQuA::gse50660_equa
  )
  expected_rows <- c(gse40279_equa = 656L, gse50660_equa = 464L)

  for (name in names(examples)) {
    dat <- examples[[name]]
    expect_s3_class(dat, "data.frame")
    expect_equal(nrow(dat), expected_rows[[name]])
    expect_named(
      dat,
      c("sample_id", "chronological_age", "estimated_age")
    )
    expect_false(anyNA(dat))
    expect_equal(anyDuplicated(dat$sample_id), 0L)
    expect_true(all(grepl("^GSM[0-9]+$", dat$sample_id)))
    expect_true(all(is.finite(dat$chronological_age)))
    expect_true(all(is.finite(dat$estimated_age)))
  }

  expect_equal(range(examples$gse40279_equa$chronological_age), c(19, 101))
  expect_equal(range(examples$gse50660_equa$chronological_age), c(38, 67))
})

test_that("public GEO example can run through the package interface", {
  dat <- EQuA::gse40279_equa
  calibration <- dat[seq_len(500), ]
  new_data <- dat[501:506, ]

  fit <- fit_equa(
    chronological_age = calibration$chronological_age,
    estimated_age = calibration$estimated_age,
    method = "kernel_end",
    probs = c(0.1, 0.5, 0.9)
  )
  profile <- predict(
    fit,
    new_estimated_age = new_data$estimated_age,
    chronological_age = new_data$chronological_age
  )

  expect_s3_class(profile, "equa_profile")
  expect_equal(dim(profile$qea), c(6L, 3L))
  expect_true(all(is.finite(profile$qea)))
  expect_true(all(is.finite(profile$qeaa)))
})

test_that("installed GEO provenance manifest is available", {
  manifest_file <- system.file(
    "extdata",
    "gse_example_manifest.csv",
    package = "EQuA"
  )
  expect_true(nzchar(manifest_file))
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)

  expect_setequal(
    manifest$object_name,
    c("gse40279_equa", "gse50660_equa")
  )
  expect_equal(manifest$package_n, c(656L, 464L))
  expect_true(all(nzchar(manifest$upstream_sha256)))
})
