# Public GEO example data quality review

## Intended use and grain

`gse40279_equa` and `gse50660_equa` are example inputs for EQuA. Each row is
one public GEO sample accession. The retained interface is deliberately
limited to:

- `sample_id`: public GSM accession;
- `chronological_age`: age in years from GEO sample metadata;
- `estimated_age`: upstream Horvath v1 DNA methylation age in years.

They are not raw-data objects, clock-calculation datasets, or reproductions of
the original publications' analyses.

## Validation results

| Check | GSE40279 | GSE50660 |
|:--|--:|--:|
| Rows | 656 | 464 |
| Unique GSM accessions | 656 | 464 |
| Duplicate rows | 0 | 0 |
| Missing retained values | 0 | 0 |
| Chronological-age range (years) | 19–101 | 38–67 |
| Estimated-age range (years) | 16.083–95.224 | 36.629–76.210 |
| CA/EA Pearson correlation | 0.9183 | 0.7656 |
| Horvath model-CpG coverage | 353/353 (100%) | 351/353 (99.43%) |

All sample identifiers match `^GSM[0-9]+$`. Row counts and age ranges agree
with the corresponding GEO records and publications. Retained values are
finite.

## Lineage checks

The package generation script requires the following upstream tables and
fails when their MD5 checksum changes:

| Accession | Frozen analysis-ready SHA-256 | MD5 |
|:--|:--|:--|
| GSE40279 | `0aa46d15bcdd8ef7cb2de644298e9736bb21d2accdc1337a8c7db56bcb3cf175` | `cc9b746f43176aede437619a20eb0526` |
| GSE50660 | `4689d047720c7ad4f843edf2c516f60baa027e4dde6170f11bd2ce6d1dd36c72` | `5adf18afcc3a6166d258b44ba31454fa` |

The frozen upstream validation report independently records 656 and 464 rows,
unique sample keys, matching schemas, and no numerical differences after
rebuilding the two source tables.

## Known limitations

- The distributed objects contain derived Horvath v1 estimates, not the full
  methylation matrices. They cannot independently reproduce clock scoring.
- GSE50660 lacks two of the 353 Horvath model CpGs. Its estimates use the
  available CpGs only and should be interpreted with that limitation.
- The raw GEO series-matrix checksum was not preserved in the compact package
  lineage. The checksum-controlled parent for the package objects is therefore
  the frozen analysis-ready table, while stable GEO URLs document the raw
  source.
- No claim is made that these datasets are representative calibration
  populations for a new scientific analysis. They are included as public
  software examples.
