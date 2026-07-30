#' Fit an EQuA calibration model
#'
#' Fits a post-hoc EQuA calibration model using chronological age, an existing
#' scalar estimated age, and optional numeric conditioning variables.
#'
#' @param chronological_age Numeric vector of chronological ages, in the same
#'   units as `estimated_age`.
#' @param estimated_age Numeric vector of scalar ages from an existing
#'   estimator, in the same units as `chronological_age`.
#' @param covariates Optional numeric vector, matrix, or data frame. Required
#'   only for `method = "kernel_ex"`.
#' @param method One of `"kernel_end"` (End-EQuA), `"po"` (Po-EQuA), or
#'   `"kernel_ex"` (Ex-EQuA).
#' @param probs Strictly increasing percentile levels between zero and one.
#' @param bandwidth Optional bandwidth vector in the original units of the
#'   conditioning variables. The default applies the documented rule separately
#'   to every conditioning dimension.
#' @param na_action Either `"omit"` or `"fail"`.
#' @param support_action Default prediction behavior in low-support regions:
#'   `"warn"`, `"na"`, `"error"`, or the explicitly requested `"nearest"`
#'   fallback.
#' @param min_effective_n Positive threshold used to flag low effective kernel
#'   sample size.
#' @param nearest_k Number of neighbors used only when
#'   `support_action = "nearest"`.
#'
#' @return An object of class `equa_fit`. Kernel fits retain the
#'   calibration-level quantities needed for later prediction.
#'
#' @section Data privacy:
#' An `equa_fit` object retains calibration-level EAD and conditioning values
#' required for prediction. It may therefore contain participant-level
#' information even though identifiers are not stored.
#'
#' @examples
#' calibration <- subset(
#'   simulated_equa_data,
#'   sample_role == "calibration"
#' )
#' fit <- fit_equa(
#'   chronological_age = calibration$chronological_age,
#'   estimated_age = calibration$estimated_age,
#'   method = "kernel_end",
#'   probs = c(0.25, 0.5, 0.75)
#' )
#' fit
#'
#' @export
fit_equa <- function(
    chronological_age,
    estimated_age,
    covariates = NULL,
    method = c("kernel_end", "po", "kernel_ex"),
    probs = seq(0.1, 0.9, by = 0.1),
    bandwidth = NULL,
    na_action = c("omit", "fail"),
    support_action = c("warn", "na", "error", "nearest"),
    min_effective_n = 10,
    nearest_k = 30) {
  method <- .match_equa_method(method)
  probs <- .validate_probs(probs)
  na_action <- match.arg(na_action)
  support <- .validate_support_settings(
    support_action,
    min_effective_n,
    nearest_k
  )
  inputs <- .prepare_calibration_inputs(
    chronological_age,
    estimated_age,
    covariates,
    method,
    na_action
  )

  ead <- inputs$chronological_age - inputs$estimated_age
  conditioning <- switch(
    method,
    po = matrix(
      numeric(),
      nrow = length(ead),
      ncol = 0L,
      dimnames = list(NULL, character())
    ),
    kernel_end = matrix(
      inputs$estimated_age,
      ncol = 1L,
      dimnames = list(NULL, "estimated_age")
    ),
    kernel_ex = cbind(
      estimated_age = inputs$estimated_age,
      inputs$covariates
    )
  )

  transform <- .conditioning_transform(conditioning)
  bandwidth_info <- .resolve_bandwidth(bandwidth, transform)

  object <- list(
    method = method,
    probs = probs,
    n_input = inputs$n_input,
    n_calibration = length(ead),
    n_omitted = inputs$n_omitted,
    calibration_ead = ead,
    calibration_conditioning = transform$scaled,
    conditioning_names = colnames(conditioning),
    conditioning_center = transform$center,
    conditioning_scale = transform$scale,
    conditioning_range = transform$range,
    bandwidth = bandwidth_info$raw,
    bandwidth_scaled = bandwidth_info$scaled,
    bandwidth_selected = bandwidth_info$selected,
    na_action = na_action,
    support_action = support$support_action,
    min_effective_n = support$min_effective_n,
    nearest_k = support$nearest_k,
    call = match.call(),
    privacy_notice = paste(
      "This object retains calibration-level values required for prediction;",
      "review disclosure risk before saving or sharing it."
    )
  )
  class(object) <- "equa_fit"
  object
}
