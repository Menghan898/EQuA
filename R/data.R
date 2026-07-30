#' Simulated age-estimator data for EQuA examples
#'
#' A small synthetic dataset reserved for package documentation, interface
#' examples, and non-sensitive tests. It does not contain participant data and
#' is not intended to reproduce the JSM 2026 analyses.
#'
#' @format A data frame with 240 rows and 5 variables:
#' \describe{
#'   \item{sample_id}{Synthetic sample identifier.}
#'   \item{sample_role}{Either `"calibration"` or `"prediction"`.}
#'   \item{chronological_age}{Simulated chronological age in years.}
#'   \item{estimated_age}{Simulated scalar estimated age in years.}
#'   \item{z}{Simulated continuous conditioning variable.}
#' }
#'
#' @source Generated deterministically by
#' `data-raw/generate_simulated_equa_data.R`.
"simulated_equa_data"

#' GSE40279 inputs for EQuA examples
#'
#' Compact, analysis-ready inputs derived from the public GSE40279 whole-blood
#' DNA methylation study. The object contains all 656 public GEO samples and is
#' intended for examples of fitting and cross-fitting EQuA profiles. It does
#' not contain methylation beta values.
#'
#' `estimated_age` is the Horvath v1 multi-tissue DNA methylation age in years,
#' calculated upstream using coefficients exported from `methylclockData`
#' 1.14.0 and a local implementation matching `methylclock::anti.trafo`. All
#' 353 model CpGs were available; no imputation was used. These derived scores
#' are not age estimates reported by the original GSE40279 publication.
#'
#' @format A data frame with 656 rows and 3 variables:
#' \describe{
#'   \item{sample_id}{Public GEO sample accession (GSM).}
#'   \item{chronological_age}{Chronological age in years.}
#'   \item{estimated_age}{Horvath v1 estimated age in years.}
#' }
#'
#' @source GEO series
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279};
#' Hannum et al. (2013), \doi{10.1016/j.molcel.2012.10.016};
#' Horvath (2013), \doi{10.1186/gb-2013-14-10-r115}. Full derivation metadata
#' are installed in `extdata/gse_example_manifest.csv` and documented in
#' `development/data-provenance.md`.
"gse40279_equa"

#' GSE50660 inputs for EQuA examples
#'
#' Compact, analysis-ready inputs derived from the public GSE50660 peripheral-
#' blood DNA methylation study. The object contains all 464 public GEO samples
#' and is intended for examples of fitting and cross-fitting EQuA profiles. It
#' does not contain methylation beta values.
#'
#' `estimated_age` is the Horvath v1 multi-tissue DNA methylation age in years,
#' calculated upstream using coefficients exported from `methylclockData`
#' 1.14.0 and a local implementation matching `methylclock::anti.trafo`. Of
#' 353 model CpGs, 351 were available (99.43 percent). The score used the
#' available CpGs only; cohort-mean and K-nearest-neighbor imputation were not
#' used. These derived scores are not age estimates reported by the original
#' GSE50660 publication.
#'
#' @format A data frame with 464 rows and 3 variables:
#' \describe{
#'   \item{sample_id}{Public GEO sample accession (GSM).}
#'   \item{chronological_age}{Chronological age in years.}
#'   \item{estimated_age}{Horvath v1 estimated age in years.}
#' }
#'
#' @source GEO series
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50660};
#' Tsaprouni et al. (2014), \doi{10.4161/15592294.2014.969637};
#' Horvath (2013), \doi{10.1186/gb-2013-14-10-r115}. Full derivation metadata
#' are installed in `extdata/gse_example_manifest.csv` and documented in
#' `development/data-provenance.md`.
"gse50660_equa"
