# Generate the non-sensitive package example data.
#
# Run this script from the package root. The generated object is documentation
# scaffolding only and does not encode an EQuA estimator.

set.seed(20260730)

n <- 240L
chronological_age <- stats::runif(n, min = 20, max = 85)
z <- stats::rnorm(n)
age_dependent_bias <- 0.04 * (chronological_age - 52.5)
estimated_age <- chronological_age - age_dependent_bias +
  1.5 * z + stats::rnorm(n, sd = 4)

simulated_equa_data <- data.frame(
  sample_id = sprintf("sim_%03d", seq_len(n)),
  sample_role = rep(c("calibration", "prediction"), c(200L, 40L)),
  chronological_age = round(chronological_age, 3),
  estimated_age = round(estimated_age, 3),
  z = round(z, 4),
  stringsAsFactors = FALSE
)

dir.create("data", showWarnings = FALSE)
save(
  simulated_equa_data,
  file = "data/simulated_equa_data.rda",
  compress = "xz"
)
