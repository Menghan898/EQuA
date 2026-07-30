test_that("weighted quantiles use the left-continuous empirical inverse", {
  result <- EQuA:::.weighted_quantile_left(
    x = c(1, 2, 3),
    weights = c(0.2, 0.3, 0.5),
    probs = c(0.1, 0.2, 0.2000001, 0.5, 0.5000001, 0.99)
  )

  expect_equal(result, c(1, 1, 2, 2, 3, 3))
})

test_that("weighted quantiles preserve ties and omit unusable weights", {
  result <- EQuA:::.weighted_quantile_left(
    x = c(2, 1, 1, NA, 100),
    weights = c(0.7, 0.1, 0.2, 1, 0),
    probs = c(0.25, 0.3, 0.3000001, 0.9)
  )

  expect_equal(result, c(1, 1, 2, 2))
  expect_true(all(is.na(
    EQuA:::.weighted_quantile_left(1:3, c(0, NA, -1), c(0.2, 0.8))
  )))
})

test_that("weighted quantiles match the frozen reference definition", {
  set.seed(18)
  x <- sample(c(-2:5, NA_real_), 80, replace = TRUE)
  w <- c(runif(75), 0, -1, NA, Inf, 1e300)
  probs <- seq(0.01, 0.99, length.out = 41)

  expect_equal(
    EQuA:::.weighted_quantile_left(x, w, probs),
    reference_weighted_quantile(x, w, probs)
  )
})
