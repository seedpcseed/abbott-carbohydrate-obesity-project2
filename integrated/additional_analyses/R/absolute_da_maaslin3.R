#!/usr/bin/env Rscript
# qPCR-scaled genus DA with MaAsLin 3 (BioC 3.23 / R >= 4.6).

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else Sys.getenv(
  "ABBOTT_ROOT",
  unset = "/home/patcseed/projects/abbott-carbohydrate-obesity-project2"
)
lib <- file.path(root, "integrated/additional_analyses/renv_bioc323/library")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

out_dir <- file.path(root, "integrated/results/additional_analyses_2026-07-19/absolute_da")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("maaslin3", quietly = TRUE)) {
  writeLines("maaslin3 is not installed", file.path(out_dir, "maaslin3_BLOCKED.md"))
  quit(status = 1)
}
suppressPackageStartupMessages(library(maaslin3))

meta <- read_csv(
  file.path(root, "integrated/results/additional_analyses_2026-07-19/tables/genus_da_sample_meta.csv"),
  show_col_types = FALSE
)
rel <- read_csv(
  file.path(root, "integrated/results/additional_analyses_2026-07-19/tables/genus_relative_for_da.csv"),
  show_col_types = FALSE
)

stopifnot(identical(meta$sample_key, rel$sample_key))

joined <- meta %>%
  filter(!is.na(gene_copies_per_ul), gene_copies_per_ul > 0) %>%
  mutate(
    obesity_group = factor(
      ifelse(as.character(group) == "case", "obesity", "healthy_weight"),
      levels = c("healthy_weight", "obesity")
    ),
    carbohydrate = factor(carbohydrate),
    timepoint = factor(timepoint),
    donor_id = factor(as.character(donor_id)),
    aliquot_id = factor(as.character(aliquot_id)),
    read_depth = as.numeric(asv_read_depth)
  ) %>%
  inner_join(rel, by = "sample_key")

message("MaAsLin3 samples: ", nrow(joined), " donors: ", dplyr::n_distinct(joined$donor_id))

genus_cols <- setdiff(names(rel), "sample_key")
prev_ok <- vapply(genus_cols, function(g) {
  tmp <- joined %>%
    transmute(donor_id, carbohydrate, timepoint, present = .data[[g]] > 0) %>%
    group_by(carbohydrate, timepoint) %>%
    summarise(
      donor_prev = n_distinct(donor_id[present]) / n_distinct(donor_id),
      .groups = "drop"
    )
  any(tmp$donor_prev >= 0.25, na.rm = TRUE)
}, logical(1))
keep_genera <- genus_cols[prev_ok]
message("Genera passing donor prevalence filter: ", length(keep_genera), " / ", length(genus_cols))

feature_mat <- as.matrix(joined[, keep_genera, drop = FALSE])
rownames(feature_mat) <- joined$sample_key
scale_vec <- joined$gene_copies_per_ul
feature_scaled <- sweep(feature_mat, 1, scale_vec, `*`)

verify_rows <- character()
fusi_cols <- grep("Fusicatenibacter", colnames(feature_scaled), ignore.case = TRUE, value = TRUE)
if (length(fusi_cols)) {
  fcol <- fusi_cols[[1]]
  i <- which(joined$carbohydrate == "SDC" & joined$timepoint %in% c("48H", "48"))[1]
  if (!is.na(i)) {
    verify_rows <- c(
      sprintf("sample=%s", joined$sample_key[i]),
      sprintf("rel=%g", feature_mat[i, fcol]),
      sprintf("gene_copies=%g", scale_vec[i]),
      sprintf("scaled=%g", feature_scaled[i, fcol]),
      sprintf("product_check=%g", feature_mat[i, fcol] * scale_vec[i])
    )
  }
}
writeLines(verify_rows, file.path(out_dir, "maaslin3_qpcr_scale_verify.txt"))

metadata <- as.data.frame(joined %>%
  select(sample_key, carbohydrate, timepoint, obesity_group, donor_id, read_depth))
rownames(metadata) <- metadata$sample_key
metadata$sample_key <- NULL

# Sanitize genus names that break model assembly (e.g. asterisks).
colnames(feature_scaled) <- make.names(colnames(feature_scaled), unique = TRUE)

tmpdir <- file.path(out_dir, "maaslin3_fit")
unlink(tmpdir, recursive = TRUE)
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

# Donor random intercept only: donor+aliquot random effects triggered a
# post-fit assembly error ("supplied row-names must match list extent") in
# maaslin3 1.4.0 on this table. Aliquot remains in metadata for diagnostics.
fit <- tryCatch({
  maaslin3(
    input_data = feature_scaled,
    input_metadata = metadata,
    output = tmpdir,
    formula = "~ carbohydrate * timepoint + obesity_group + (1|donor_id)",
    normalization = "NONE",
    transform = "LOG",
    median_comparison_abundance = FALSE,
    warn_prevalence = FALSE,
    max_significance = 0.25,
    cores = 1L,
    plot_summary_plot = FALSE,
    plot_associations = FALSE
  )
}, error = function(e) {
  writeLines(conditionMessage(e), file.path(out_dir, "maaslin3_ERROR.txt"))
  e
})

res_files <- list.files(tmpdir, pattern = "all_results|significant", recursive = TRUE, full.names = TRUE)
file.copy(res_files, out_dir, overwrite = TRUE)
if (!inherits(fit, "error") && !is.null(fit$fit_data_abundance$results)) {
  write_csv(as.data.frame(fit$fit_data_abundance$results), file.path(out_dir, "maaslin3_abundance_results.csv"))
}
if (!inherits(fit, "error") && !is.null(fit$fit_data_prevalence$results)) {
  write_csv(as.data.frame(fit$fit_data_prevalence$results), file.path(out_dir, "maaslin3_prevalence_results.csv"))
}
write_csv(tibble(genus = keep_genera), file.path(out_dir, "maaslin3_genera_kept.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "maaslin3_sessionInfo.txt"))
unlink(file.path(out_dir, "maaslin3_BLOCKED.md"), force = TRUE)
unlink(file.path(out_dir, "BLOCKED.md"), force = TRUE)
message("MaAsLin3 finished; outputs in ", out_dir)
