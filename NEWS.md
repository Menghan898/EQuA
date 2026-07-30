# EQuA 0.0.0.9000

## Initial development implementation

- Added the initial package skeleton, documentation layout, citation metadata,
  and continuous-integration configuration.
- Implemented Po-EQuA, End-EQuA, and Ex-EQuA through
  `fit_equa()` and `predict.equa_fit()`.
- Implemented repeated chronological-age-stratified K-fold profiles through
  `crossfit_equa()`.
- Added S3 printing, summaries, plotting, data-frame conversion, and explicit
  support diagnostics.
- Fixed the reference definitions for left-continuous weighted quantiles,
  dimension-dependent bandwidth selection, fold construction, missing values,
  unit-equivalent scaling, and repeated-profile aggregation.
- Added regression tests against independent calculations derived from the
  internal reproducibility programs.
- Added compact full-cohort examples `gse40279_equa` and `gse50660_equa`,
  together with checksum-controlled generation, installed machine-readable
  provenance, original-study and Horvath-clock citations, explicit
  missing-CpG policies, and a data quality review.
- Replaced the interim all-rights-reserved notice with the MIT License.
- Added a compact two-dataset GEO clock-QC table to the repository README.
- Renamed the simulated-data role `"application"` to `"prediction"` and use
  `new_data` consistently for prediction examples.
- Confirmed Menghan Yi as the sole software author and maintainer for the
  initial development release; poster authorship remains separately recorded.
