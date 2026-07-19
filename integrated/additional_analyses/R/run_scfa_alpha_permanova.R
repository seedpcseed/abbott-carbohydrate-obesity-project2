#!/usr/bin/env Rscript
# Export adjusted SCFA change contrasts, donor-aware alpha models,
# and term-level PERMANOVA + betadisper.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(vegan)
  library(phyloseq)
})

root <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = FALSE)
# Prefer project root regardless of invocation cwd
args_root <- Sys.getenv("ABBOTT_ROOT", unset = "")
if (nzchar(args_root)) {
  root <- args_root
} else {
  # Walk up from this script
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
  )
  root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
}

out <- file.path(root, "integrated", "results", "additional_analyses_2026-07-19")
dir.create(file.path(out, "scfa_contrasts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "alpha_permanova"), recursive = TRUE, showWarnings = FALSE)

message("ROOT=", root)

## ---- SCFA nested models + carbohydrate change contrasts -----------------
scfa_path <- file.path(
  root, "scfa_metabolomics", "data",
  "removed_qcs_quant_results_20250722_PFBBr_PS2720_20250731.csv"
)
meta_canon <- readr::read_csv(
  file.path(root, "integrated", "metadata", "canonical_experimental_units_scfa.csv"),
  show_col_types = FALSE
)

scfa_raw <- readr::read_csv(scfa_path, show_col_types = FALSE) %>%
  dplyr::rename(sampleid = sampleid) %>%
  dplyr::filter(!stringr::str_detect(sampleid, regex("QC|Blank|Standard|pool", ignore_case = TRUE)))

scfa_long <- scfa_raw %>%
  tidyr::pivot_longer(
    cols = c(acetate, propionate, butyrate),
    names_to = "analyte",
    values_to = "concentration"
  ) %>%
  dplyr::mutate(sample_id = stringr::str_trim(tolower(sampleid))) %>%
  dplyr::inner_join(
    meta_canon %>%
      dplyr::mutate(sample_id = stringr::str_trim(tolower(sample_id))),
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    donor_id = factor(donor_id),
    aliquot_id = factor(aliquot_id),
    well_id = factor(well_id),
    group = factor(group, levels = c("control", "case")),
    carbohydrate = factor(
      carbohydrate,
      levels = c("No Added Carb", "RDC", "SDC")
    ),
    timepoint = factor(timepoint, levels = c("0H", "48H"))
  )

contrast_rows <- list()
vc_rows <- list()
for (an in c("acetate", "propionate", "butyrate")) {
  d <- dplyr::filter(scfa_long, analyte == an)
  fit <- tryCatch(
    lmer(
      concentration ~ group * carbohydrate * timepoint
        + (1 | donor_id) + (1 | aliquot_id) + (1 | well_id),
      data = d,
      REML = TRUE
    ),
    error = function(e) {
      message("Fallback model for ", an, ": ", e$message)
      lmer(
        concentration ~ group * carbohydrate * timepoint
          + (1 | donor_id) + (1 | aliquot_id),
        data = d,
        REML = TRUE
      )
    }
  )
  vc <- as.data.frame(VarCorr(fit))
  vc$analyte <- an
  vc_rows[[an]] <- vc

  # Marginalize over group; pairwise change contrasts among carbohydrates
  emm <- emmeans(fit, ~ timepoint | carbohydrate)
  # Force equal weights over group by regrid if needed: use emmeans with group averaged
  emm_g <- emmeans(fit, ~ timepoint * carbohydrate)
  # Change per carbohydrate (48 - 0), then pairwise among carbs
  # Use contrast on interaction grid
  cont <- contrast(
    emmeans(fit, ~ timepoint * carbohydrate),
    interaction = c("pairwise", "pairwise"),
    adjust = "none"
  )
  # Better: compute delta per carb then pairwise
  emm2 <- emmeans(fit, ~ timepoint | carbohydrate)
  deltas <- contrast(emm2, method = "revpairwise") # 48H - 0H within carb
  deltas_df <- as.data.frame(summary(deltas, infer = TRUE))
  deltas_df$analyte <- an
  deltas_df$contrast_family <- "within_carb_time_change"

  # Pairwise carbohydrate differences in change: use interaction contrast
  # (carb i 48 - carb i 0) - (carb j 48 - carb j 0)
  emm3 <- emmeans(fit, ~ carbohydrate * timepoint)
  # Custom contrasts for change differences
  lev <- levels(interaction(d$carbohydrate, d$timepoint, sep = ":", drop = TRUE))
  # Build through emmeans pairs of deltas using contrast combinations
  # Use emmeans::contrast with interaction on carbohydrate for the 48-0 contrast
  emm_delta <- contrast(emm3, interaction = "pairwise", by = NULL)
  # Simpler reliable approach: estimate change per carb, then pairwise among those
  # Extract EMMs for each carb*time
  emm_tab <- as.data.frame(summary(emm3, infer = TRUE))
  get_mean <- function(carb, tp) {
    emm_tab$emmean[emm_tab$carbohydrate == carb & emm_tab$timepoint == tp][1]
  }
  get_var_row <- function(carb, tp) {
    # Use vcov of emmeans
    NULL
  }
  # Use contrast objects from emmeans for DiD among carbs
  # Create a custom contrast matrix via pairs on change
  change_contrasts <- list(
    "RDC_change" = c(
      "No Added Carb:0H" = 0, "No Added Carb:48H" = 0,
      "RDC:0H" = -1, "RDC:48H" = 1,
      "SDC:0H" = 0, "SDC:48H" = 0
    ),
    "SDC_change" = c(
      "No Added Carb:0H" = 0, "No Added Carb:48H" = 0,
      "RDC:0H" = 0, "RDC:48H" = 0,
      "SDC:0H" = -1, "SDC:48H" = 1
    ),
    "NC_change" = c(
      "No Added Carb:0H" = -1, "No Added Carb:48H" = 1,
      "RDC:0H" = 0, "RDC:48H" = 0,
      "SDC:0H" = 0, "SDC:48H" = 0
    ),
    "RDC_vs_NC_change" = c(
      "No Added Carb:0H" = 1, "No Added Carb:48H" = -1,
      "RDC:0H" = -1, "RDC:48H" = 1,
      "SDC:0H" = 0, "SDC:48H" = 0
    ),
    "SDC_vs_NC_change" = c(
      "No Added Carb:0H" = 1, "No Added Carb:48H" = -1,
      "RDC:0H" = 0, "RDC:48H" = 0,
      "SDC:0H" = -1, "SDC:48H" = 1
    ),
    "SDC_vs_RDC_change" = c(
      "No Added Carb:0H" = 0, "No Added Carb:48H" = 0,
      "RDC:0H" = 1, "RDC:48H" = -1,
      "SDC:0H" = -1, "SDC:48H" = 1
    )
  )
  # Align names with emmeans grid
  grid <- as.data.frame(emm3@grid)
  grid$key <- paste(grid$carbohydrate, grid$timepoint, sep = ":")
  built <- lapply(change_contrasts, function(cc) {
    v <- setNames(rep(0, nrow(grid)), grid$key)
    v[names(cc)] <- cc
    as.numeric(v)
  })
  names(built) <- names(change_contrasts)
  ct <- contrast(emm3, method = built, adjust = "tukey")
  ct_df <- as.data.frame(summary(ct, infer = TRUE))
  ct_df$analyte <- an
  ct_df$contrast_family <- "carb_change_and_pairwise"
  # Within-family Tukey already applied to the 6 contrasts; also BH within analyte pairwise subset
  contrast_rows[[an]] <- dplyr::bind_rows(deltas_df, ct_df)
}

contrasts_all <- dplyr::bind_rows(contrast_rows)
readr::write_csv(contrasts_all, file.path(out, "scfa_contrasts", "carbohydrate_change_contrasts.csv"))
readr::write_csv(dplyr::bind_rows(vc_rows), file.path(out, "scfa_contrasts", "scfa_variance_components.csv"))
message("Wrote SCFA contrasts")

## ---- Alpha diversity donor-aware models ---------------------------------
# Build from genus-relative? Prefer phyloseq if available in cache; else estimate from ASV
asv <- readr::read_csv(
  file.path(
    root, "microbiome_analysis", "zr24558.16S_250813.zymo",
    "00...AllSamples.Bac16Sv34", "DADA2_ASV_Distribution", "ASV_Abundance_Table.csv"
  ),
  show_col_types = FALSE
)
sample_col <- names(asv)[1]
mat <- as.matrix(asv[, -1, drop = FALSE])
rownames(mat) <- asv[[sample_col]]
# Diversity metrics
shannon <- vegan::diversity(mat, index = "shannon")
observed <- rowSums(mat > 0)
simpson <- vegan::diversity(mat, index = "simpson")
alpha <- tibble::tibble(
  customer_label = names(shannon),
  Shannon = as.numeric(shannon),
  Observed = as.numeric(observed),
  Simpson = as.numeric(simpson)
) %>%
  dplyr::mutate(
    customer_label = stringr::str_replace_all(stringr::str_trim(customer_label), "\\.+$", ""),
    customer_label = stringr::str_replace_all(customer_label, "\\.\\.", ".")
  )

parse_lab <- function(lab) {
  m <- stringr::str_match(lab, "^(\\d+[A-Za-z])\\.(R[12]|S[12]|NC)\\.*(0H|48H)$")
  if (any(is.na(m))) return(tibble::tibble())
  aliq <- toupper(m[, 2])
  cond <- toupper(m[, 3])
  tp <- toupper(m[, 4])
  carb <- dplyr::case_when(
    stringr::str_starts(cond, "R") ~ "RDC",
    stringr::str_starts(cond, "S") ~ "SDC",
    TRUE ~ "No Added Carb"
  )
  well <- ifelse(cond == "NC", 1L, as.integer(stringr::str_sub(cond, 2, 2)))
  donor <- stringr::str_extract(aliq, "^\\d+")
  tibble::tibble(
    customer_label = lab,
    donor_id = donor,
    aliquot_id = aliq,
    carbohydrate = carb,
    well_repeat = well,
    timepoint = tp
  )
}
lab_meta <- dplyr::bind_rows(lapply(alpha$customer_label, parse_lab))
alpha_m <- alpha %>%
  dplyr::inner_join(lab_meta, by = "customer_label") %>%
  dplyr::mutate(donor_id = as.character(donor_id)) %>%
  dplyr::left_join(
    meta_canon %>%
      dplyr::mutate(donor_id = as.character(donor_id)) %>%
      dplyr::distinct(donor_id, group),
    by = "donor_id"
  ) %>%
  dplyr::filter(!is.na(group))

# Aggregate to donor × carb × time
alpha_donor <- alpha_m %>%
  dplyr::group_by(donor_id, group, carbohydrate, timepoint) %>%
  dplyr::summarise(
    Shannon = mean(Shannon, na.rm = TRUE),
    Observed = mean(Observed, na.rm = TRUE),
    Simpson = mean(Simpson, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    donor_id = factor(donor_id),
    group = factor(group, levels = c("control", "case")),
    carbohydrate = factor(carbohydrate, levels = c("No Added Carb", "RDC", "SDC")),
    timepoint = factor(timepoint, levels = c("0H", "48H"))
  )

alpha_model_rows <- list()
for (metric in c("Shannon", "Observed", "Simpson")) {
  form <- as.formula(paste0(metric, " ~ group * carbohydrate * timepoint + (1 | donor_id)"))
  fit <- lmer(form, data = alpha_donor, REML = TRUE)
  an <- as.data.frame(anova(fit, type = 3))
  an$metric <- metric
  an$term <- rownames(an)
  alpha_model_rows[[metric]] <- an
  # emmeans for group at 48H No Added Carb
}
alpha_anova <- dplyr::bind_rows(alpha_model_rows)
readr::write_csv(alpha_anova, file.path(out, "alpha_permanova", "alpha_mixed_anova.csv"))
readr::write_csv(alpha_donor, file.path(out, "alpha_permanova", "alpha_donor_cell.csv"))

# Prespecified contrasts: Healthy vs Obese at 48H under each carb
alpha_contrast_rows <- list()
for (metric in c("Shannon", "Observed", "Simpson")) {
  fit <- lmer(
    as.formula(paste0(metric, " ~ group * carbohydrate * timepoint + (1 | donor_id)")),
    data = alpha_donor, REML = TRUE
  )
  emm <- emmeans(fit, ~ group | carbohydrate * timepoint)
  ct <- contrast(emm, method = "revpairwise", adjust = "BH")
  df <- as.data.frame(summary(ct, infer = TRUE))
  df$metric <- metric
  alpha_contrast_rows[[metric]] <- df
}
readr::write_csv(
  dplyr::bind_rows(alpha_contrast_rows),
  file.path(out, "alpha_permanova", "alpha_group_contrasts.csv")
)
message("Wrote alpha models")

## ---- PERMANOVA term-level + betadisper ----------------------------------
# Aggregate relative genus profile? Use ASV Bray-Curtis on donor×carb×time means of relative ASV
rel <- sweep(mat, 1, rowSums(mat), "/")
rel[!is.finite(rel)] <- 0
rel_df <- as.data.frame(rel)
rel_df$customer_label <- rownames(rel)
rel_m <- rel_df %>%
  dplyr::mutate(
    customer_label = stringr::str_replace_all(stringr::str_trim(customer_label), "\\.+$", ""),
    customer_label = stringr::str_replace_all(customer_label, "\\.\\.", ".")
  ) %>%
  dplyr::inner_join(lab_meta, by = "customer_label") %>%
  dplyr::mutate(donor_id = as.character(donor_id)) %>%
  dplyr::left_join(
    meta_canon %>%
      dplyr::mutate(donor_id = as.character(donor_id)) %>%
      dplyr::distinct(donor_id, group),
    by = "donor_id"
  ) %>%
  dplyr::filter(!is.na(group), !is.na(carbohydrate), !is.na(timepoint))

taxa_cols <- setdiff(names(rel_m), c("customer_label", names(lab_meta), "group"))
donor_comm <- rel_m %>%
  dplyr::group_by(donor_id, group, carbohydrate, timepoint) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(taxa_cols), mean), .groups = "drop")

comm_mat <- as.matrix(donor_comm[, taxa_cols])
rownames(comm_mat) <- paste(donor_comm$donor_id, donor_comm$carbohydrate, donor_comm$timepoint, sep = "::")
meta_comm <- donor_comm %>%
  dplyr::select(donor_id, group, carbohydrate, timepoint) %>%
  dplyr::mutate(
    group = factor(group),
    carbohydrate = factor(carbohydrate),
    timepoint = factor(timepoint),
    donor_id = factor(donor_id)
  )

dist_all <- vegan::vegdist(comm_mat, method = "bray")
set.seed(20260719)
perm_all <- vegan::adonis2(
  dist_all ~ group * carbohydrate + timepoint,
  data = meta_comm,
  permutations = 9999,
  by = "terms"
)
perm_all_df <- as.data.frame(perm_all)
perm_all_df$term <- rownames(perm_all_df)
perm_all_df$model <- "all_timepoints_donor_aggregated"
readr::write_csv(perm_all_df, file.path(out, "alpha_permanova", "permanova_terms_all.csv"))

# 48H only
idx48 <- meta_comm$timepoint == "48H"
dist48 <- vegan::vegdist(comm_mat[idx48, , drop = FALSE], method = "bray")
meta48 <- meta_comm[idx48, , drop = FALSE]
perm48 <- vegan::adonis2(
  dist48 ~ group * carbohydrate,
  data = meta48,
  permutations = 9999,
  by = "terms"
)
perm48_df <- as.data.frame(perm48)
perm48_df$term <- rownames(perm48_df)
perm48_df$model <- "48H_donor_aggregated"
readr::write_csv(perm48_df, file.path(out, "alpha_permanova", "permanova_terms_48h.csv"))

# betadisper for carbohydrate and group at 48H
disp_rows <- list()
for (fac in c("carbohydrate", "group")) {
  bd <- vegan::betadisper(dist48, meta48[[fac]])
  pt <- vegan::permutest(bd, permutations = 9999)
  disp_rows[[fac]] <- tibble::tibble(
    factor = fac,
    F = as.numeric(pt$tab$F[1]),
    pvalue = as.numeric(pt$tab$`Pr(>F)`[1]),
    n = nrow(meta48)
  )
}
readr::write_csv(dplyr::bind_rows(disp_rows), file.path(out, "alpha_permanova", "betadisper_48h.csv"))
message("Wrote PERMANOVA + dispersion")

## ---- Manifest ------------------------------------------------------------
manifest <- tibble::tibble(
  file = list.files(out, recursive = TRUE),
  path = file.path(out, list.files(out, recursive = TRUE))
) %>%
  dplyr::mutate(size_bytes = file.info(path)$size)
readr::write_csv(manifest, file.path(out, "manifest.csv"))
message("Done")
