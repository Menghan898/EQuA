# EQuA data provenance

## Principles

Every dataset distributed with EQuA must have:

- a clear analytical purpose;
- a documented source and retrieval date;
- reproducible preparation code;
- explicit inclusion and exclusion rules;
- variable definitions and units;
- licensing and attribution review;
- no restricted participant information.

Package data contain inputs for examples, not precomputed qEA or qEAA results.

## Simulated example data

`simulated_equa_data` is generated deterministically by
`data-raw/generate_simulated_equa_data.R`.

It contains only synthetic identifiers, chronological age, estimated age, a
continuous conditioning variable, and a calibration/prediction role. It does
not reproduce or summarize any real cohort.

## Public GEO examples

The package includes two compact data frames:

- `gse40279_equa`: 656 samples from GSE40279;
- `gse50660_equa`: 464 samples from GSE50660.

Each row represents one public GEO sample accession, and each object retains
only:

- `sample_id`: GSM accession;
- `chronological_age`: chronological age in years;
- `estimated_age`: an upstream Horvath v1 multi-tissue DNA methylation age in
  years.

The objects do not include methylation beta matrices, original study-specific
identifiers, sex, site, plate, ethnicity, smoking status, other clock scores,
precomputed EAD, qEA, or qEAA. EQuA does not calculate the retained source
clock.

### GSE40279

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279>
- Series matrix:
  <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz>
- Original study: Hannum et al. (2013), *Genome-wide methylation profiles
  reveal quantitative views of human aging rates*, Molecular Cell 49(2):
  359–367.
- PMID: 23177740
- DOI: <https://doi.org/10.1016/j.molcel.2012.10.016>
- Original and retained sample count: 656/656
- Chronological-age range: 19–101 years
- Horvath model-CpG coverage: 353/353; no imputation

### GSE50660

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50660>
- Series matrix:
  <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50660/matrix/GSE50660_series_matrix.txt.gz>
- Original study: Tsaprouni et al. (2014), *Cigarette smoking reduces DNA
  methylation levels at multiple genomic loci but the effect is partially
  reversible upon cessation*.
- PMID: 25424692
- DOI: <https://doi.org/10.4161/15592294.2014.969637>
- Original and retained sample count: 464/464
- Chronological-age range: 38–67 years
- Horvath model-CpG coverage: 351/353 (99.43%); available-CpG-only scoring,
  without cohort-mean or K-nearest-neighbor imputation

### Estimated-age derivation

The `estimated_age` column is a derived Horvath v1 score and was not reported
as such by either original study. The upstream workflow:

1. read public beta values and chronological-age metadata from the GEO series
   matrix;
2. used Horvath v1 coefficients exported from `methylclockData` 1.14.0;
3. formed the linear score from available model CpGs;
4. applied a local implementation matching `methylclock::anti.trafo`;
5. joined the score to GEO metadata by GSM accession;
6. built checksum-validated analysis-ready clock tables; and
7. retained only GSM, chronological age, and Horvath v1 estimated age for the
   package objects.

The Horvath clock source is Horvath (2013), *DNA methylation age of human
tissues and cell types*, Genome Biology 14:R115, PMID 24138928,
<https://doi.org/10.1186/gb-2013-14-10-r115>.

### Rebuild and integrity controls

`data-raw/generate_geo_example_data.R` requires `EQUA_REPRO_ROOT` to point to
the frozen upstream project. It refuses to generate package data unless each
analysis-ready parent table has the expected row count, schema, sample-key
properties, and MD5 checksum. The parent SHA-256 checksums are:

- GSE40279:
  `0aa46d15bcdd8ef7cb2de644298e9736bb21d2accdc1337a8c7db56bcb3cf175`;
- GSE50660:
  `4689d047720c7ad4f843edf2c516f60baa027e4dde6170f11bd2ce6d1dd36c72`.

The upstream GEO matrices were inspected on 2026-07-01. Their raw compressed
checksums were not retained in the compact package lineage, so the frozen
analysis-ready tables—not a newly downloaded GEO file—are the
checksum-controlled direct parents of the package data. This limitation is
recorded rather than silently substituting a current GEO checksum.

The installed machine-readable record is available at:

```r
system.file(
  "extdata",
  "gse_example_manifest.csv",
  package = "EQuA"
)
```

Detailed range, uniqueness, missingness, and lineage checks are recorded in
`development/geo-example-data-quality.md`.

### Attribution and redistribution boundary

The GEO records are public and remain attributable to their original
investigators and publications. EQuA redistributes only compact derived
CA/EA example inputs, with stable accession links and citations. It does not
redistribute the full GEO matrices or the Horvath coefficient table. Inclusion
in EQuA does not replace any citation or reuse requirements attached to the
original studies, GEO, or the clock implementation.
