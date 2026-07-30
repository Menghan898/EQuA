# EQuA development plan

## Status

This repository now contains a tested pure-R reference implementation.
Po-EQuA, Kernel End-EQuA, Kernel Ex-EQuA, prediction diagnostics, and repeated
cross-fitting are implemented. The public interface and governance metadata
remain development-stage.

## A. Scope and non-goals

The first development release accepts chronological age, an existing
estimated age, and optional numeric conditioning variables. It returns qEA
and qEAA profiles using Po-EQuA, kernel End-EQuA, or kernel Ex-EQuA.

The package will not calculate source age estimators, preprocess methylation
data, or contain cohort-specific analyses.

## B. Statistical specification

The following conventions are frozen for the current reference version:

1. left-continuous weighted empirical quantiles without interpolation;
2. the reproduction-program bandwidth rule, selected inside each training
   fold;
3. chronological-age-stratified fold construction;
4. method-specific complete cases;
5. pointwise median aggregation across repeats by default;
6. effective sample size below 10 as the default low-support flag;
7. errors for zero-variance conditioning variables;
8. calibration-SD Euclidean distance for the nearest-distance diagnostic;
9. unit-equivalent internal scaling;
10. no implicit nearest-neighbor fallback.

The exact equations and edge cases are maintained in
`vignettes/method-specification.Rmd`.

## C. Public API and object model

The implemented public interface is deliberately small:

- `fit_equa()`;
- `predict.equa_fit()`;
- `crossfit_equa()`;
- `print()`, `summary()`, `plot()`, and `as.data.frame()` methods.

The implemented S3 classes are:

- `equa_fit`;
- `equa_profile`;
- `equa_crossfit`, inheriting from `equa_profile`.

Cross-fitting defaults to `keep_repeats = FALSE`. Fold assignments and
fold-level diagnostics will be configurable separately.

## D. Implementation principles

- Begin with a transparent pure-R reference implementation.
- Keep dependencies light; core calculations should rely on base R and
  `stats`.
- Do not export low-level numerical helpers.
- Add compiled code only after the reference implementation and tests are
  stable.
- Preserve results under pure changes of conditioning-variable units.
- Never silently switch to a different estimator in low-support regions.
- Avoid retaining subject identifiers or unused participant variables.

## E. Validation strategy

Before public analytical use, tests will cover:

- hand-verifiable weighted quantiles;
- monotonicity over percentile levels;
- unit invariance;
- deterministic cross-fitting under a fixed seed;
- strict out-of-fold prediction;
- missing and non-finite inputs;
- zero-variance conditioning variables;
- support diagnostics and action semantics;
- privacy-safe printing, summaries, errors, and warnings;
- agreement between any optimized implementation and the pure-R reference.

README synchronization is tested separately in CI by rebuilding `README.md`
from `README.Rmd` and requiring a clean Git diff.

## F. Data and release roadmap

1. **Completed locally — skeleton:** metadata, citations, documentation, CI,
   and simulated data.
2. **Completed locally — reference implementation:** Po-, End-, and Ex-EQuA
   fit and prediction.
3. **Completed locally — cross-fitting and diagnostics:** repeated out-of-fold
   profiles and support reporting.
4. **Completed locally — public examples:** curated GSE40279 and GSE50660
   inputs with checksum-controlled generation, installed provenance, and data
   quality review.
5. **Performance:** optional compiled optimization after numerical equivalence
   is established.

## Governance decisions required before public release

- confirm the package name and capitalization;
- approve the software and poster citation forms;

Menghan Yi is the confirmed maintainer, with public contact email
`menghany@umich.edu`. Menghan Yi selected the MIT License for the software on
2026-07-30 and confirmed that the initial development release should list
Menghan Yi as its sole software author (`aut`) and maintainer (`cre`). The JSM
2026 poster retains its separate six-author citation. Additional software
contributors can be added later when their package-development roles are
established.
