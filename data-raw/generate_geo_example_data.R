# Generate compact public GEO inputs for package examples.
#
# This script intentionally does not download methylation matrices or calculate
# epigenetic clocks. It extracts a validated CA/EA interface from the frozen
# upstream analysis-ready tables. Set EQUA_REPRO_ROOT to the root of that
# reproducibility project before running this script.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  "data-raw/generate_geo_example_data.R"
}
package_root <- normalizePath(
  file.path(dirname(script_file), ".."),
  mustWork = TRUE
)

repro_root <- Sys.getenv("EQUA_REPRO_ROOT", unset = "")
if (!nzchar(repro_root)) {
  stop(
    "Set EQUA_REPRO_ROOT to the frozen reproducibility-project root.",
    call. = FALSE
  )
}
repro_root <- normalizePath(repro_root, mustWork = TRUE)

specification <- data.frame(
  accession = c("GSE40279", "GSE50660"),
  object_name = c("gse40279_equa", "gse50660_equa"),
  expected_rows = c(656L, 464L),
  expected_md5 = c(
    "cc9b746f43176aede437619a20eb0526",
    "5adf18afcc3a6166d258b44ba31454fa"
  ),
  stringsAsFactors = FALSE
)

build_one <- function(accession, object_name, expected_rows, expected_md5) {
  source_file <- file.path(
    repro_root,
    "03_analysis_ready",
    "data",
    accession,
    "base_clock_table.csv"
  )
  if (!file.exists(source_file)) {
    stop("Missing upstream table: ", source_file, call. = FALSE)
  }

  observed_md5 <- unname(tools::md5sum(source_file))
  if (!identical(observed_md5, expected_md5)) {
    stop(
      "Checksum mismatch for ", accession, ": expected ", expected_md5,
      ", observed ", observed_md5, ".",
      call. = FALSE
    )
  }

  upstream <- utils::read.csv(
    source_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required <- c(
    "dataset",
    "sample_id",
    "chronological_age_years",
    "horvathv1"
  )
  missing_columns <- setdiff(required, names(upstream))
  if (length(missing_columns)) {
    stop(
      "Missing required columns for ", accession, ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(upstream) != expected_rows) {
    stop(
      "Unexpected row count for ", accession, ": expected ", expected_rows,
      ", observed ", nrow(upstream), ".",
      call. = FALSE
    )
  }
  if (!all(upstream$dataset == accession)) {
    stop("Unexpected dataset label in ", accession, ".", call. = FALSE)
  }

  result <- data.frame(
    sample_id = as.character(upstream$sample_id),
    chronological_age = as.numeric(upstream$chronological_age_years),
    estimated_age = as.numeric(upstream$horvathv1),
    stringsAsFactors = FALSE
  )

  if (anyNA(result) || any(!is.finite(result$chronological_age)) ||
      any(!is.finite(result$estimated_age))) {
    stop("Missing or non-finite retained values in ", accession, ".", call. = FALSE)
  }
  if (anyDuplicated(result$sample_id)) {
    stop("Duplicate GSM accessions in ", accession, ".", call. = FALSE)
  }
  if (!all(grepl("^GSM[0-9]+$", result$sample_id))) {
    stop("Invalid GSM accession in ", accession, ".", call. = FALSE)
  }

  assign(object_name, result, envir = .GlobalEnv)
  output_file <- file.path(package_root, "data", paste0(object_name, ".rda"))
  save(
    list = object_name,
    file = output_file,
    envir = .GlobalEnv,
    compress = "xz",
    version = 2
  )
  message("Wrote ", output_file)
  invisible(result)
}

for (i in seq_len(nrow(specification))) {
  build_one(
    accession = specification$accession[[i]],
    object_name = specification$object_name[[i]],
    expected_rows = specification$expected_rows[[i]],
    expected_md5 = specification$expected_md5[[i]]
  )
}
