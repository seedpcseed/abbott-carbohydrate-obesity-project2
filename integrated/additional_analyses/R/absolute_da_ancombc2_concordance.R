#!/usr/bin/env Rscript
# Relative-abundance ANCOM-BC2 coefficient export for concordance with MaAsLin3/ALDEx3.
# Uses the render environment (ANCOMBC 2.4.x); input is matched-well genus counts.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else Sys.getenv(
  "ABBOTT_ROOT",
  unset = "/home/patcseed/projects/abbott-carbohydrate-obesity-project2"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(phyloseq)
  library(ANCOMBC)
})

out_dir <- file.path(root, "integrated/results/additional_analyses_2026-07-19/absolute_da")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260719)

counts <- read_csv(
  file.path(root, "integrated/results/additional_analyses_2026-07-19/tables/genus_counts_matched.csv"),
  show_col_types = FALSE
)
meta <- read_csv(
  file.path(root, "integrated/results/additional_analyses_2026-07-19/tables/genus_da_sample_meta.csv"),
  show_col_types = FALSE
)
stopifnot(identical(counts$sample_key, meta$sample_key))

feature_cols <- setdiff(names(counts), "sample_key")
otu <- t(as.matrix(counts[, feature_cols, drop = FALSE]))
colnames(otu) <- counts$sample_key
# Drop empty taxa
otu <- otu[rowSums(otu) > 0, , drop = FALSE]

sd <- meta %>%
  transmute(
    sample_key,
    carbohydrate = factor(carbohydrate),
    timepoint = factor(timepoint),
    obesity_group = factor(
      ifelse(group == "case", "obesity", "healthy_weight"),
      levels = c("healthy_weight", "obesity")
    ),
    donor_id = factor(as.character(donor_id)),
    # ANCOMBC 2.4 in this render env does not parse `*` / `:` interaction
    # tokens in fix_formula; encode carb×time as an explicit factor.
    carb_time = factor(paste(carbohydrate, timepoint, sep = "_"))
  ) %>%
  as.data.frame()
rownames(sd) <- sd$sample_key
sd$sample_key <- NULL
sd <- sd[colnames(otu), , drop = FALSE]

tax <- data.frame(Genus = rownames(otu), row.names = rownames(otu))
ps <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  sample_data(sd),
  tax_table(as.matrix(tax))
)

message("ANCOM-BC2 samples=", nsamples(ps), " genera=", ntaxa(ps), " donors=", nlevels(sd$donor_id))

fit <- tryCatch({
  ancombc2(
    data = ps,
    tax_level = "Genus",
    fix_formula = "carb_time + obesity_group",
    rand_formula = "(1 | donor_id)",
    p_adj_method = "BH",
    prv_cut = 0.10,
    lib_cut = 1000,
    s0_perc = 0.05,
    group = "carb_time",
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = 1,
    verbose = TRUE
  )
}, error = function(e) e)

if (inherits(fit, "error")) {
  writeLines(conditionMessage(fit), file.path(out_dir, "ancombc2_ERROR.txt"))
  # Fallback without random effect
  fit <- tryCatch({
    ancombc2(
      data = ps,
      tax_level = "Genus",
      fix_formula = "carb_time + obesity_group",
      rand_formula = NULL,
      p_adj_method = "BH",
      prv_cut = 0.10,
      lib_cut = 1000,
      s0_perc = 0.05,
      group = "carb_time",
      struc_zero = FALSE,
      neg_lb = FALSE,
      alpha = 0.05,
      n_cl = 1,
      verbose = TRUE
    )
  }, error = function(e) e)
}

if (inherits(fit, "error")) {
  writeLines(conditionMessage(fit), file.path(out_dir, "ancombc2_ERROR.txt"))
  quit(status = 1)
}

res <- as.data.frame(fit$res)
write_csv(res, file.path(out_dir, "ancombc2_relative_results_wide.csv"))

# Long form for concordance: one row per taxon × effect
lfc_cols <- grep("^lfc_", names(res), value = TRUE)
se_cols <- grep("^se_", names(res), value = TRUE)
q_cols <- grep("^q_", names(res), value = TRUE)
taxon_col <- if ("taxon" %in% names(res)) "taxon" else names(res)[1]

long <- lapply(lfc_cols, function(lc) {
  effect <- sub("^lfc_", "", lc)
  se_c <- paste0("se_", effect)
  q_c <- paste0("q_", effect)
  tibble(
    feature = as.character(res[[taxon_col]]),
    effect = effect,
    ancom_lfc = res[[lc]],
    ancom_se = if (se_c %in% names(res)) res[[se_c]] else NA_real_,
    ancom_q = if (q_c %in% names(res)) res[[q_c]] else NA_real_
  )
}) %>% bind_rows()

write_csv(long, file.path(out_dir, "ancombc2_relative_results_long.csv"))

# Three-method concordance for shared SDC×48H / obesity terms where names align
m <- read_csv(file.path(out_dir, "maaslin3_abundance_results.csv"), show_col_types = FALSE)
a <- read_csv(file.path(out_dir, "aldex3_results_svar_0.25.csv"), show_col_types = FALSE)

map_maaslin <- m %>%
  transmute(
    feature = as.character(feature),
    effect = dplyr::case_when(
      grepl("SDC:timepoint48H|SDC:timepoint 48H", paste(metadata, value), ignore.case = TRUE) ~ "carbohydrateSDC:timepoint48H",
      grepl("RDC:timepoint48H", paste(metadata, value), ignore.case = TRUE) ~ "carbohydrateRDC:timepoint48H",
      metadata == "timepoint" & value == "48H" ~ "timepoint48H",
      metadata == "obesity_group" ~ "obesity_groupobesity",
      metadata == "carbohydrate" & value == "SDC" ~ "carbohydrateSDC",
      metadata == "carbohydrate" & value == "RDC" ~ "carbohydrateRDC",
      TRUE ~ paste0(metadata, value)
    ),
    maaslin_coef = coef,
    maaslin_q = qval_individual,
    maaslin_joint_q = qval_joint
  )

map_aldex <- a %>%
  transmute(
    feature = as.character(entity),
    effect = as.character(parameter),
    aldex_est = estimate,
    aldex_q = p.val.adj
  )

map_ancom <- long %>%
  mutate(
    effect = dplyr::case_when(
      grepl("SDC_48H|SDC.*48H", effect, ignore.case = TRUE) ~ "carbohydrateSDC:timepoint48H",
      grepl("RDC_48H|RDC.*48H", effect, ignore.case = TRUE) ~ "carbohydrateRDC:timepoint48H",
      grepl("No Added Carb_48H|No.Added.Carb_48H|No_Added_Carb_48H", effect, ignore.case = TRUE) ~ "timepoint48H",
      grepl("obesity", effect, ignore.case = TRUE) ~ "obesity_groupobesity",
      grepl("SDC_0H|SDC.*0H", effect, ignore.case = TRUE) ~ "carbohydrateSDC",
      grepl("RDC_0H|RDC.*0H", effect, ignore.case = TRUE) ~ "carbohydrateRDC",
      TRUE ~ effect
    )
  )

conc <- map_maaslin %>%
  full_join(map_aldex, by = c("feature", "effect")) %>%
  full_join(map_ancom, by = c("feature", "effect")) %>%
  mutate(
    dir_maaslin = sign(maaslin_coef),
    dir_aldex = sign(aldex_est),
    dir_ancom = sign(ancom_lfc),
    concordant_maaslin_aldex = !is.na(dir_maaslin) & !is.na(dir_aldex) & dir_maaslin == dir_aldex,
    concordant_all3 = !is.na(dir_maaslin) & !is.na(dir_aldex) & !is.na(dir_ancom) &
      dir_maaslin == dir_aldex & dir_maaslin == dir_ancom
  )

write_csv(conc, file.path(out_dir, "concordance_ancombc2_maaslin3_aldex3.csv"))

# Focused primary interaction summary
focus <- conc %>%
  filter(effect %in% c("carbohydrateSDC:timepoint48H", "carbohydrateRDC:timepoint48H", "obesity_groupobesity", "timepoint48H"))
write_csv(focus, file.path(out_dir, "concordance_focus_effects.csv"))

fusi <- focus %>% filter(grepl("Fusi", feature, ignore.case = TRUE))
write_csv(fusi, file.path(out_dir, "concordance_fusicatenibacter.csv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "ancombc2_sessionInfo.txt"))
message("ANCOM-BC2 + concordance written to ", out_dir)
message("Focus Fusicatenibacter rows: ", nrow(fusi))
