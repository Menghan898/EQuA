# Data preparation

Files in this directory generate documented package datasets. The directory is
excluded from the built source package but remains version controlled.

Current scripts:

- `generate_simulated_equa_data.R`: creates the deterministic, non-sensitive
  `simulated_equa_data` object.
- `generate_geo_example_data.R`: extracts the compact `gse40279_equa` and
  `gse50660_equa` objects from checksum-verified upstream analysis-ready
  tables. Set `EQUA_REPRO_ROOT` to the frozen upstream project root before
  running it.

The installed, machine-readable provenance record for the public GEO examples
is `inst/extdata/gse_example_manifest.csv`. It records source accessions,
stable URLs, original publications, retained variables, clock scoring details,
feature coverage, missing-CpG policies, and upstream-table checksums.
