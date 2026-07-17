#' Build a Table 1-style clinical characteristics table (HLCA format)
#'
#' Collapses cell-level metadata from a Seurat object into a per-sample
#' clinical characteristics table, strictly following the column structure
#' of the Human Lung Cell Atlas (HLCA) reference Table 1. Uses an explicit
#' whitelist mechanism to avoid accidentally concatenating cell-level
#' variables (e.g. \code{percent.mt}, \code{nCount_RNA}, \code{Phase},
#' \code{seurat_clusters}) into meaningless per-sample strings.
#'
#' @param seurat_obj A \code{Seurat} object whose \code{meta.data} contains
#'   one row per cell, with a sample identifier column (default \code{"sample"}).
#' @param sample_col Character. Name of the column in \code{meta.data} that
#'   identifies the sample/library each cell belongs to. Default \code{"sample"}.
#' @param sample_level_cols Character vector. Whitelist of column names in
#'   \code{meta.data} that are true sample/subject-level variables (i.e.
#'   constant within a sample) and should therefore be collapsed to one
#'   value per sample. Any column not listed here is treated as cell-level
#'   and excluded from the collapse step. Defaults to the standard set used
#'   in the LungAgingERV project (see Details).
#' @param audit_suffix Character. Suffix used to identify "locked/audited"
#'   columns (e.g. \code{age_excel_final_v3}) that should be used to fill
#'   in gaps in the primary columns. Default \code{"_v3"}.
#' @param subject_id_parser A function that takes a single sample name
#'   (character scalar) and returns the parsed \code{subject_ID} (character
#'   scalar, \code{NA_character_} if unparseable). Defaults to
#'   \code{\link{parse_subject_id_default}}, which handles the
#'   \code{Hs_WT_*} and \code{Hs_Banovich_GSE*} naming conventions used in
#'   LungAgingERV. Supply a custom function for other projects
#'   (HeartERVmap, PH-ERVAtlas, etc.) with different sample naming schemes.
#' @param hlca_core_pattern Character (regex). Samples whose
#'   \code{dataset_origin} matches this pattern are labelled
#'   \code{"core"} in the output \code{HLCA_core_or_extension} column;
#'   all others are labelled \code{"project_new (not part of HLCA)"}.
#'   Default \code{"Banovich"}.
#' @param save_output Logical. If \code{TRUE}, writes the result table to
#'   \code{file.path(output_dir, paste0(output_prefix, ".csv"))} and
#'   \code{.xlsx} (with QC sheets). Default \code{FALSE} (function returns
#'   results in-memory only; safer default for a package function).
#' @param output_dir Character. Directory to write outputs to, if
#'   \code{save_output = TRUE}. Default \code{"."}.
#' @param output_prefix Character. File name prefix (without extension) for
#'   saved outputs. Default \code{"TableS1_HLCA_format"}.
#'
#' @return A list with components:
#' \describe{
#'   \item{table1}{The final per-sample table, arranged and renamed to
#'     match the HLCA reference Table 1 column structure.}
#'   \item{sample_summary}{The intermediate per-sample collapsed table
#'     (whitelist columns only, before HLCA-format renaming), useful for
#'     downstream merging (e.g. with \code{merge_table1_sources()}).}
#'   \item{na_report}{Data frame summarising missingness per column in
#'     \code{table1}.}
#'   \item{excluded_cols}{Character vector of cell-level columns that were
#'     explicitly excluded from the collapse step.}
#'   \item{n_cells_check}{A named list with \code{table1_total} and
#'     \code{seurat_total}, and a logical \code{match} flag, from the
#'     internal QC check that \code{n_cells} sums to the total number of
#'     cells in \code{seurat_obj}.}
#' }
#'
#' @details
#' The default \code{sample_level_cols} whitelist is:
#' \code{dataset_origin, AgeGroup, AgeGroup_short, age, Biological_Sex,
#' single_cell_platform, sequencing_platform, cell_ranger_version,
#' tissue_dissociation_protocol, harmonized_ethnicity,
#' selfreported_ethnicity_as_collected, smoking_status, BMI,
#' cause_of_death, fresh_or_frozen, anatomical_region_level_1/2/3,
#' cells_or_nuclei, tissue_sampling_type}, plus the \code{_v3}-suffixed
#' audit columns (\code{age_excel_final_v3}, \code{sex_final_v3},
#' \code{AgeGroup_final_v3}, \code{Self_Reported_Ethnicity_v3},
#' \code{Ever_Smoker_v3}, \code{Disease_Status_v3}).
#'
#' If any whitelisted column is found to have more than one distinct
#' non-missing value within a single sample, a warning is raised (the
#' first value is kept) — this usually indicates the column was
#' mistakenly included in the whitelist and is actually cell-level.
#'
#' @examples
#' \dontrun{
#' res <- build_table1_hlca_format(immune.combined)
#' res$table1
#' res$na_report
#'
#' # Custom subject_ID parser for a different project naming scheme:
#' my_parser <- function(s) sub("^Sample_", "", s)
#' res2 <- build_table1_hlca_format(immune.combined, subject_id_parser = my_parser)
#' }
#'
#' @export
build_table1_hlca_format <- function(
    seurat_obj,
    sample_col = "sample",
    sample_level_cols = c(
      "dataset_origin", "AgeGroup", "AgeGroup_short", "age", "Biological_Sex",
      "single_cell_platform", "sequencing_platform", "cell_ranger_version",
      "tissue_dissociation_protocol", "harmonized_ethnicity",
      "selfreported_ethnicity_as_collected", "smoking_status", "BMI",
      "cause_of_death", "fresh_or_frozen",
      "anatomical_region_level_1", "anatomical_region_level_2", "anatomical_region_level_3",
      "cells_or_nuclei", "tissue_sampling_type",
      "age_excel_final_v3", "sex_final_v3", "AgeGroup_final_v3",
      "Self_Reported_Ethnicity_v3", "Ever_Smoker_v3", "Disease_Status_v3"
    ),
    audit_suffix = "_v3",
    subject_id_parser = parse_subject_id_default,
    hlca_core_pattern = "Banovich",
    save_output = FALSE,
    output_dir = ".",
    output_prefix = "TableS1_HLCA_format"
) {

  stopifnot(inherits(seurat_obj, "Seurat"))
  meta <- seurat_obj@meta.data
  if (!sample_col %in% colnames(meta)) {
    stop(sprintf("sample_col '%s' not found in seurat_obj@meta.data", sample_col))
  }
  if (sample_col != "sample") {
    meta[["sample"]] <- meta[[sample_col]]
  }

  # ----------------------------------------------------------------
  # Step 1. Whitelist detection (cell-level columns are explicitly excluded)
  # ----------------------------------------------------------------
  sample_level_cols_available <- intersect(sample_level_cols, colnames(meta))
  missing_from_whitelist <- setdiff(sample_level_cols, colnames(meta))
  if (length(missing_from_whitelist) > 0) {
    message("The following whitelisted columns were not found in meta.data and will be skipped: ",
            paste(missing_from_whitelist, collapse = ", "))
  }

  cell_level_cols_excluded <- setdiff(colnames(meta), c("sample", sample_level_cols_available))

  # ----------------------------------------------------------------
  # Step 2. Collapse: keep the unique non-NA value per sample per column
  # ----------------------------------------------------------------
  safe_first_unique <- function(x, colname) {
    x <- x[!is.na(x) & as.character(x) != ""]
    u <- unique(as.character(x))
    if (length(u) == 0) return(NA_character_)
    if (length(u) > 1) {
      warning(sprintf(
        "Column '%s' has %d distinct non-missing values within a single sample: %s (keeping the first; check whether this column truly belongs in the whitelist)",
        colname, length(u), paste(u, collapse = " | ")
      ))
      return(u[1])
    }
    u
  }

  sample_summary <- meta %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sample_level_cols_available),
                     ~ safe_first_unique(.x, dplyr::cur_column())),
      n_cells_this_project = dplyr::n(),
      .groups = "drop"
    )

  # ----------------------------------------------------------------
  # Step 3. Parse subject_ID via injectable parser function
  # ----------------------------------------------------------------
  sample_summary$subject_ID <- vapply(sample_summary$sample, subject_id_parser,
                                       FUN.VALUE = character(1))

  # ----------------------------------------------------------------
  # Step 4. Field harmonization: primary column first, audit(_v3) column as fallback
  # ----------------------------------------------------------------
  audit_col <- function(base) paste0(base, audit_suffix)
  has_col <- function(nm) nm %in% colnames(sample_summary)

  sample_summary <- sample_summary %>%
    dplyr::mutate(
      age_final = dplyr::coalesce(
        suppressWarnings(as.numeric(.data[["age"]])),
        if (has_col(audit_col("age_excel_final"))) suppressWarnings(as.numeric(.data[[audit_col("age_excel_final")]])) else NA_real_
      ),
      sex_final = tolower(dplyr::coalesce(
        .data[["Biological_Sex"]],
        if (has_col(audit_col("sex_final"))) .data[[audit_col("sex_final")]] else NA_character_
      )),
      ethnicity_final = dplyr::coalesce(
        .data[["harmonized_ethnicity"]],
        if (has_col(audit_col("Self_Reported_Ethnicity"))) .data[[audit_col("Self_Reported_Ethnicity")]] else NA_character_
      ),
      selfreported_ethnicity_final = dplyr::coalesce(
        .data[["selfreported_ethnicity_as_collected"]],
        if (has_col(audit_col("Self_Reported_Ethnicity"))) .data[[audit_col("Self_Reported_Ethnicity")]] else NA_character_
      ),
      smoking_status_final = dplyr::coalesce(
        .data[["smoking_status"]],
        if (has_col(audit_col("Ever_Smoker"))) {
          dplyr::case_when(
            .data[[audit_col("Ever_Smoker")]] == "Y" ~ "active_or_former",
            .data[[audit_col("Ever_Smoker")]] == "N" ~ "never",
            TRUE ~ NA_character_
          )
        } else NA_character_
      ),
      lung_condition_final = "Healthy"  # default assumption for healthy-aging cohorts; adjust for disease projects
    )

  # ----------------------------------------------------------------
  # Step 5. HLCA_core_or_extension + dataset/study labels
  # ----------------------------------------------------------------
  sample_summary <- sample_summary %>%
    dplyr::mutate(
      HLCA_core_or_extension = dplyr::if_else(
        stringr::str_detect(.data[["dataset_origin"]], hlca_core_pattern),
        "core", "project_new (not part of HLCA)"
      ),
      dataset = .data[["dataset_origin"]],
      study   = .data[["dataset_origin"]],
      mixed_harmonized_ethnicity = NA_character_,
      mixed_selfreported_ethnicity_as_collected = NA_character_
    )

  # ----------------------------------------------------------------
  # Step 6. Reassemble in HLCA reference column order
  # ----------------------------------------------------------------
  table1 <- sample_summary %>%
    dplyr::transmute(
      sample                    = .data$sample,
      HLCA_core_or_extension    = .data$HLCA_core_or_extension,
      dataset                   = .data$dataset,
      study                     = .data$study,
      n_cells                   = .data$n_cells_this_project,
      subject_ID                = .data$subject_ID,
      `age*`                    = .data$age_final,
      sex                       = .data$sex_final,
      lung_condition            = .data$lung_condition_final,
      cells_or_nuclei           = .data[["cells_or_nuclei"]],
      single_cell_platform      = .data[["single_cell_platform"]],
      tissue_sampling_type      = .data[["tissue_sampling_type"]],
      tissue_dissociation_protocol = .data[["tissue_dissociation_protocol"]],
      harmonized_ethnicity      = .data$ethnicity_final,
      selfreported_ethnicity_as_collected = .data$selfreported_ethnicity_final,
      mixed_harmonized_ethnicity = .data$mixed_harmonized_ethnicity,
      mixed_selfreported_ethnicity_as_collected = .data$mixed_selfreported_ethnicity_as_collected,
      smoking_status            = .data$smoking_status_final,
      BMI                       = suppressWarnings(as.numeric(.data[["BMI"]])),
      cause_of_death            = .data[["cause_of_death"]],
      sequencing_platform       = .data[["sequencing_platform"]],
      cell_ranger_version       = .data[["cell_ranger_version"]],
      fresh_or_frozen           = .data[["fresh_or_frozen"]],
      anatomical_region_level_1 = .data[["anatomical_region_level_1"]],
      anatomical_region_level_2 = .data[["anatomical_region_level_2"]],
      anatomical_region_level_3 = .data[["anatomical_region_level_3"]]
    ) %>%
    dplyr::arrange(.data$HLCA_core_or_extension, .data$`age*`, .data$sample)

  # ----------------------------------------------------------------
  # Step 7. QC: n_cells sum vs actual total cell count
  # ----------------------------------------------------------------
  total_table1 <- sum(table1$n_cells)
  total_seurat <- ncol(seurat_obj)
  n_cells_check <- list(
    table1_total  = total_table1,
    seurat_total  = total_seurat,
    match         = identical(total_table1, total_seurat)
  )
  if (!n_cells_check$match) {
    warning(sprintf(
      "n_cells sum in table1 (%d) does not match total cells in seurat_obj (%d). Check for unsampled cells.",
      total_table1, total_seurat
    ))
  }

  na_report <- data.frame(
    column = colnames(table1),
    n_missing = sapply(table1, function(x) sum(is.na(x))),
    row.names = NULL
  )
  na_report$pct_missing <- round(100 * na_report$n_missing / nrow(table1), 1)
  na_report <- na_report[order(-na_report$n_missing), ]

  result <- list(
    table1          = table1,
    sample_summary  = sample_summary,
    na_report       = na_report,
    excluded_cols   = cell_level_cols_excluded,
    n_cells_check   = n_cells_check
  )

  # ----------------------------------------------------------------
  # Step 8. Optional file output
  # ----------------------------------------------------------------
  if (isTRUE(save_output)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    csv_path  <- file.path(output_dir, paste0(output_prefix, ".csv"))
    xlsx_path <- file.path(output_dir, paste0(output_prefix, ".xlsx"))

    utils::write.csv(table1, csv_path, row.names = FALSE)

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "TableS1_HLCA_format")
    openxlsx::writeData(wb, "TableS1_HLCA_format", table1)
    openxlsx::addWorksheet(wb, "Missing_value_report")
    openxlsx::writeData(wb, "Missing_value_report", na_report)
    openxlsx::addWorksheet(wb, "Excluded_cell_level_cols")
    openxlsx::writeData(wb, "Excluded_cell_level_cols",
                         data.frame(excluded_column = cell_level_cols_excluded))
    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

    message("Saved: ", csv_path, " and ", xlsx_path)
  }

  result
}


#' Default subject_ID parser for LungAgingERV sample naming conventions
#'
#' Parses \code{subject_ID} from sample names following either the
#' \code{Hs_WT_{AgeGroup}_{subjectID}_{age}} pattern (CPI_WT cohort) or
#' the \code{Hs_Banovich_{GSExxxxx}_{AgeGroup}_{subjectID}} pattern
#' (Banovich/HLCA-extension cohort).
#'
#' @param s Character scalar. A single sample name.
#' @return Character scalar: the parsed subject ID, or \code{NA_character_}
#'   if the sample name does not match either pattern.
#' @export
parse_subject_id_default <- function(s) {
  m_wt <- stringr::str_match(s, "^Hs_WT_(Young|Middle|Old)_([A-Za-z0-9]+)_([0-9]+)$")
  if (!is.na(m_wt[1, 1])) return(m_wt[1, 3])
  m_ban <- stringr::str_match(s, "^Hs_Banovich_(GSE[0-9]+)_(Young|Middle|Old)_([A-Za-z0-9_]+)$")
  if (!is.na(m_ban[1, 1])) return(m_ban[1, 4])
  NA_character_
}
