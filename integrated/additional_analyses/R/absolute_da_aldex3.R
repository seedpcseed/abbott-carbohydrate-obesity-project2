#!/usr/bin/env Rscript
# ALDEx3 qPCR-scale-uncertainty sensitivity (CRAN 1.2.x).
# Uses sample.sm with s.mu = log2(gene_copies_per_ul) and a provisional s.var grid.

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
set.seed(20260719)

if (!requireNamespace("ALDEx3", quietly = TRUE)) {
  writeLines("ALDEx3 not installed", file.path(out_dir, "aldex3_BLOCKED.md"))
  quit(status = 1)
}
suppressPackageStartupMessages(library(ALDEx3))

counts <- read_csv(
  file.path(root, "integrated/results/additional_analyses_2026-07-19/tables/genus_counts_matched.csv"),
  show_col_types = FALSE
)
meta <- read_csv(
  file.path(root, "integrated/results/additional_analyses_2026-07-19/tables/genus_da_sample_meta.csv"),
  show_col_types = FALSE
)
stopifnot(identical(counts$sample_key, meta$sample_key))

# features × samples (ALDEx3 convention)
feature_cols <- setdiff(names(counts), "sample_key")
Y0 <- t(as.matrix(counts[, feature_cols, drop = FALSE]))
colnames(Y0) <- counts$sample_key

# ~75% sparsity amalgamated into Other
zero_frac <- rowMeans(Y0 == 0)
keep <- names(zero_frac)[zero_frac <= 0.75]
drop <- names(zero_frac)[zero_frac > 0.75]
Y <- Y0[keep, , drop = FALSE]
if (length(drop)) {
  Y <- rbind(Y, other = colSums(Y0[drop, , drop = FALSE]))
}

X <- as.data.frame(meta %>%
  transmute(
    carbohydrate = factor(carbohydrate),
    timepoint = factor(timepoint),
    obesity_group = factor(
      ifelse(group == "case", "obesity", "healthy_weight"),
      levels = c("healthy_weight", "obesity")
    ),
    donor_id = factor(as.character(donor_id)),
    gene_copies_per_ul = gene_copies_per_ul
  ))
rownames(X) <- meta$sample_key
X <- X[colnames(Y), ]

s.mu <- log2(pmax(X$gene_copies_per_ul, 1))
names(s.mu) <- rownames(X)

# Use blmm for speed on full matrix; repeat primary s.var=0.25 with lme4
# on a top-prevalence subset if blmm succeeds.
s_var_grid <- c(0.05, 0.25, 1.0)
summaries <- list()

for (sv in s_var_grid) {
  message("ALDEx3 sample.sm s.var=", sv, " method=blmm")
  s.var <- rep(sv, length(s.mu))
  names(s.var) <- names(s.mu)
  fit <- tryCatch({
    aldex(
      Y,
      ~ carbohydrate * timepoint + obesity_group + (1 | donor_id),
      X,
      method = "blmm",
      nsample = 128,
      n.cores = 1L,
      scale = sample.sm,
      s.mu = s.mu,
      s.var = s.var
    )
  }, error = function(e) e)

  if (inherits(fit, "error")) {
    writeLines(
      c(paste0("s.var=", sv), conditionMessage(fit)),
      file.path(out_dir, sprintf("aldex3_ERROR_svar_%.2f.txt", sv))
    )
    next
  }
  tab <- as.data.frame(summary(fit))
  tab$s_var <- sv
  tab$engine <- "blmm"
  summaries[[as.character(sv)]] <- tab
  write_csv(tab, file.path(out_dir, sprintf("aldex3_results_svar_%.2f.csv", sv)))
}

# Exact lme4 sensitivity at primary s.var on a trimmed feature set for runtime.
message("ALDEx3 sample.sm s.var=0.25 method=lme4 (trimmed)")
keep_n <- min(40L, nrow(Y) - 1L)
Y_trim <- rbind(Y[seq_len(keep_n), , drop = FALSE], other = colSums(Y[setdiff(seq_len(nrow(Y)), seq_len(keep_n)), , drop = FALSE]))
s.var <- rep(0.25, length(s.mu))
names(s.var) <- names(s.mu)
fit_lme4 <- tryCatch({
  aldex(
    Y_trim,
    ~ carbohydrate * timepoint + obesity_group + (1 | donor_id),
    X,
    method = "lme4",
    nsample = 64,
    n.cores = 1L,
    scale = sample.sm,
    s.mu = s.mu,
    s.var = s.var
  )
}, error = function(e) e)
if (!inherits(fit_lme4, "error")) {
  tab <- as.data.frame(summary(fit_lme4))
  tab$s_var <- 0.25
  tab$engine <- "lme4_trimmed"
  write_csv(tab, file.path(out_dir, "aldex3_results_svar_0.25_lme4_trimmed.csv"))
  summaries[["lme4_025"]] <- tab
} else {
  writeLines(conditionMessage(fit_lme4), file.path(out_dir, "aldex3_ERROR_lme4_trimmed.txt"))
}

if (length(summaries)) {
  bind_rows(summaries) %>% write_csv(file.path(out_dir, "aldex3_results_all_svar.csv"))
}

writeLines(c(
  "Scale model: sample.sm",
  "s.mu = log2(gene_copies_per_ul)",
  "s.var grid = 0.05, 0.25, 1.0 (provisional; facility qPCR technical variance unavailable)",
  "Formula: ~ carbohydrate * timepoint + obesity_group + (1 | donor_id)",
  "Primary engine = blmm (full table); exact lme4 sensitivity on trimmed features at s.var=0.25",
  "nsample = 128 (blmm) / 64 (lme4 trimmed)",
  "Counts = round(genus_relative * asv_read_depth); genera with >75% zeros amalgamated to 'other'"
), file.path(out_dir, "aldex3_scale_notes.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "aldex3_sessionInfo.txt"))
unlink(file.path(out_dir, "aldex3_BLOCKED.md"), force = TRUE)
message("ALDEx3 finished")
