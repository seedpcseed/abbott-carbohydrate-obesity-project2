#!/usr/bin/env Rscript

# DADA2 processing pipeline (Big Data approach) for paired-end 16S FASTQ files
# - Reads raw FASTQs from `fastq_raw/`
# - Writes filtered FASTQs to `fastq_filtered/`
# - Learns error models
# - Streams per-sample dereplication, inference, and merging to per-sample seqtabs
# - Merges per-sample seqtabs in chunks, removes chimeras
# - Saves sequence tables, representative sequences, and exports an ASV table TSV
#
# This script follows the Big Data paired-end workflow:
# https://benjjneb.github.io/dada2/bigdata_paired.html
#
# Conventions in this repo:
# - Output directory: `dada2_asv/`
# - Sample columns in TSV end with `_filtered` (analysis scripts strip this suffix)

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

message_time <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

# -----------------------------
# Configuration (edit as needed)
# -----------------------------

# Directories
RAW_DIR  <- Sys.getenv("DADA2_RAW_DIR",  unset = "fastq_raw")
FILT_DIR <- Sys.getenv("DADA2_FILT_DIR", unset = "fastq_filtered")
OUT_DIR  <- Sys.getenv("DADA2_OUT_DIR",  unset = "dada2_asv")

if (!dir.exists(RAW_DIR))  stop("Missing RAW_DIR: ", RAW_DIR)
if (!dir.exists(FILT_DIR)) dir.create(FILT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR))  dir.create(OUT_DIR, recursive = TRUE)

# Optional: path to sample metadata to control sample name ordering and validation
METADATA_PATH <- Sys.getenv("DADA2_METADATA", unset = "sample.metadata.tsv")

# Filtering/trimming parameters. Adjust if you know the run specs.
# For V3V4 MiSeq runs, common defaults are below; tune as needed.
TRUNC_LEN  <- as.integer(strsplit(Sys.getenv("DADA2_TRUNC_LEN",  unset = "240,200"), ",")[[1]])
TRIM_LEFT  <- as.integer(strsplit(Sys.getenv("DADA2_TRIM_LEFT",  unset = "0,0"),   ",")[[1]])
MAX_EE     <- as.numeric(strsplit(Sys.getenv("DADA2_MAX_EE",     unset = "2,2"),    ",")[[1]])
TRUNC_Q    <- as.integer(Sys.getenv("DADA2_TRUNC_Q", unset = 2))
MAX_N      <- as.integer(Sys.getenv("DADA2_MAX_N",   unset = 0))
POOLING    <- FALSE
MULTI      <- as.logical(Sys.getenv("DADA2_MULTITHREAD", unset = "TRUE"))
VERBOSE    <- as.logical(Sys.getenv("DADA2_VERBOSE", unset = "TRUE"))

# Big Data merge chunk size for combining per-sample seqtabs
MERGE_CHUNK <- as.integer(Sys.getenv("DADA2_MERGE_CHUNK", unset = "50"))

if (length(TRUNC_LEN) != 2) stop("DADA2_TRUNC_LEN must have two comma-separated values: forward,reverse")
if (length(TRIM_LEFT) != 2) stop("DADA2_TRIM_LEFT must have two comma-separated values: forward,reverse")
if (length(MAX_EE) != 2)    stop("DADA2_MAX_EE must have two comma-separated values: forward,reverse")

# -----------------------------
# Helper: discover paired files
# -----------------------------

list_fastq_pairs <- function(raw_dir) {
  fq <- sort(list.files(raw_dir, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE))
  if (length(fq) == 0) stop("No FASTQ files found in ", raw_dir)
  fwd <- fq[grepl("_R1\\.fastq(\\.gz)?$", fq)]
  rev <- fq[grepl("_R2\\.fastq(\\.gz)?$", fq)]
  if (length(fwd) == 0 || length(rev) == 0) stop("Expected paired-end FASTQs named *_R1.fastq.gz and *_R2.fastq.gz")

  # Sample names: strip _R[12].fastq(.gz)
  get_sname <- function(x) gsub("_R[12]\\.fastq(\\.gz)?$", "", basename(x))
  sn_f <- get_sname(fwd)
  sn_r <- get_sname(rev)

  # Intersect and align
  common <- intersect(sn_f, sn_r)
  if (length(common) == 0) stop("No matching R1/R2 pairs found.")
  fwd <- fwd[match(common, sn_f)]
  rev <- rev[match(common, sn_r)]
  stopifnot(all(!is.na(fwd)), all(!is.na(rev)))

  data.frame(sample = common, fnF = fwd, fnR = rev, stringsAsFactors = FALSE)
}

pairs <- list_fastq_pairs(RAW_DIR)
message_time("Discovered ", nrow(pairs), " paired samples in ", RAW_DIR)

# Optional: if metadata exists, keep/align sample order
if (file.exists(METADATA_PATH)) {
  md <- tryCatch(readr::read_tsv(METADATA_PATH, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(md) && "sample_id" %in% names(md)) {
    wanted <- md$sample_id
    keep <- pairs$sample %in% wanted
    if (any(!keep)) {
      message_time("Dropping ", sum(!keep), " samples not present in metadata.")
      pairs <- pairs[keep, , drop = FALSE]
    }
    pairs <- pairs[match(wanted[wanted %in% pairs$sample], pairs$sample), , drop = FALSE]
  }
}
stopifnot(nrow(pairs) > 0)

# -----------------------------
# Filter and trim
# -----------------------------

filtF <- file.path(FILT_DIR, basename(pairs$fnF))
filtR <- file.path(FILT_DIR, basename(pairs$fnR))

names(filtF) <- pairs$sample
names(filtR) <- pairs$sample

message_time("Filtering and trimming reads → ", FILT_DIR)
ft <- filterAndTrim(
  fwd = pairs$fnF, filt = filtF,
  rev = pairs$fnR, filt.rev = filtR,
  truncLen = TRUNC_LEN, trimLeft = TRIM_LEFT,
  maxN = MAX_N, maxEE = MAX_EE, truncQ = TRUNC_Q,
  rm.phix = TRUE, compress = TRUE, multithread = MULTI
)

readr::write_tsv(as.data.frame(ft) %>% tibble::rownames_to_column("file"), file.path(OUT_DIR, "filter_trim_stats.tsv"))
message_time("Filtering done. Kept reads: ", sum(ft[,2]), "/", sum(ft[,1]), " (", round(100*sum(ft[,2])/sum(ft[,1]),1), "%)")

# Keep only samples with nonzero reads after filtering
keep_nonzero <- ft[,2] > 0
if (any(!keep_nonzero)) {
  message_time("Excluding ", sum(!keep_nonzero), " samples with 0 reads after filtering.")
}
pairs <- pairs[keep_nonzero, , drop = FALSE]
filtF <- filtF[keep_nonzero]
filtR <- filtR[keep_nonzero]

# -----------------------------
# Learn error rates
# -----------------------------

message_time("Learning error rates (forward)...")
errF <- learnErrors(filtF, multithread = MULTI)
message_time("Learning error rates (reverse)...")
errR <- learnErrors(filtR, multithread = MULTI)

saveRDS(errF, file.path(OUT_DIR, "errF.rds"))
saveRDS(errR, file.path(OUT_DIR, "errR.rds"))

# -----------------------------------------------
# Big Data: stream per-sample inference and merge
# -----------------------------------------------

tmp_seqtab_dir <- file.path(OUT_DIR, "per_sample_seqtab")
if (!dir.exists(tmp_seqtab_dir)) dir.create(tmp_seqtab_dir, recursive = TRUE)

samples <- pairs$sample

# Track counts per sample for read-tracking output
denoisedF_counts <- setNames(integer(length(samples)), samples)
denoisedR_counts <- setNames(integer(length(samples)), samples)
merged_counts    <- setNames(integer(length(samples)), samples)

message_time("Streaming per-sample dereplication, inference, and merging...")
for (i in seq_along(samples)) {
  s <- samples[i]
  fF <- filtF[i]
  fR <- filtR[i]

  if (!file.exists(fF) || !file.exists(fR)) {
    message_time("[WARN] Missing filtered files for ", s, "; skipping.")
    next
  }

  # Dereplicate this sample
  drF <- derepFastq(fF)
  drR <- derepFastq(fR)

  # Infer sample variants
  ddF <- dada(drF, err = errF, pool = POOLING, multithread = MULTI, verbose = VERBOSE)
  ddR <- dada(drR, err = errR, pool = POOLING, multithread = MULTI, verbose = VERBOSE)

  # Record denoised read counts
  denoisedF_counts[s] <- sum(getUniques(ddF))
  denoisedR_counts[s] <- sum(getUniques(ddR))

  # Merge pairs for this sample
  mg <- mergePairs(ddF, drF, ddR, drR, verbose = VERBOSE)
  merged_counts[s] <- if (is.null(mg) || NROW(mg) == 0) 0L else sum(mg$abundance)

  # Make per-sample sequence table and persist
  mg_list <- list(mg)
  names(mg_list) <- s
  st <- makeSequenceTable(mg_list)
  saveRDS(st, file.path(tmp_seqtab_dir, paste0(s, ".rds")))
}

# -----------------------------
# Make sequence tables and remove chimeras
# -----------------------------

message_time("Merging per-sample sequence tables in chunks of ", MERGE_CHUNK, "...")

merge_seqtab_files <- function(files, chunk_size = 50) {
  acc <- NULL
  i <- 1L
  n <- length(files)
  while (i <= n) {
    j <- min(i + chunk_size - 1L, n)
    chunk_tabs <- lapply(files[i:j], readRDS)
    chunk_merged <- Reduce(mergeSequenceTables, chunk_tabs)
    if (is.null(acc)) {
      acc <- chunk_merged
    } else {
      acc <- mergeSequenceTables(acc, chunk_merged)
    }
    i <- j + 1L
  }
  acc
}

seqtab_files <- sort(list.files(tmp_seqtab_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(seqtab_files) == 0) stop("No per-sample seqtab files were generated.")

seqtab <- merge_seqtab_files(seqtab_files, chunk_size = MERGE_CHUNK)
message_time("Sequence table dimensions: ", paste(dim(seqtab), collapse = " x "))

message_time("Removing chimeras (consensus)...")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = MULTI)
message_time("Non-chimeric table dimensions: ", paste(dim(seqtab.nochim), collapse = " x "))

saveRDS(seqtab,        file.path(OUT_DIR, "seqtab.rds"))
saveRDS(seqtab.nochim, file.path(OUT_DIR, "seqtab_nochim.rds"))

# Representative sequences (FASTA)
seqs <- colnames(seqtab)
rs_all <- Biostrings::DNAStringSet(seqs); names(rs_all) <- paste0("ASV_", seq_along(seqs))
Biostrings::writeXStringSet(rs_all, filepath = file.path(OUT_DIR, "repseqs.fasta"))

seqs_nc <- colnames(seqtab.nochim)
rs_nc <- Biostrings::DNAStringSet(seqs_nc); names(rs_nc) <- paste0("ASV_", seq_along(seqs_nc))
Biostrings::writeXStringSet(rs_nc, filepath = file.path(OUT_DIR, "repseqs.nochimera.fasta"))

# -----------------------------
# Export a tidy ASV table for downstream analysis
# -----------------------------

message_time("Exporting ASV table (TSV)...")

# Build a data frame with ASV metadata and counts
counts <- t(seqtab.nochim)  # taxa as rows

# Column naming convention: sample + "_filtered" (matches analysis scripts)
colnames(counts) <- paste0(rownames(seqtab.nochim), "_filtered")

asv_ids <- paste0("ASV_", seq_len(nrow(counts)))
seq_lengths <- nchar(rownames(counts))
total_abund <- rowSums(counts)
rel_abund   <- total_abund / sum(total_abund)

asv_table <- tibble::tibble(
  ASV_ID = asv_ids,
  Sequence = rownames(counts),
  Length = seq_lengths,
  Total_Abundance = total_abund,
  Relative_Abundance = rel_abund
) %>% dplyr::bind_cols(as.data.frame(counts))

readr::write_tsv(asv_table, file.path(OUT_DIR, "asv_table.tsv"))

# Tracking table: reads through the pipeline per sample
ft_sub <- ft[keep_nonzero, , drop = FALSE]
rownames(ft_sub) <- pairs$sample

nonchim_counts <- setNames(integer(length(samples)), samples)
nonchim_counts[rownames(seqtab.nochim)] <- rowSums(seqtab.nochim)

n_reads <- tibble::tibble(
  sample   = samples,
  input    = as.integer(ft_sub[samples, 1]),
  filtered = as.integer(ft_sub[samples, 2]),
  denoisedF = as.integer(denoisedF_counts[samples]),
  denoisedR = as.integer(denoisedR_counts[samples]),
  merged    = as.integer(merged_counts[samples]),
  nonchim   = as.integer(nonchim_counts[samples])
)
readr::write_tsv(n_reads, file.path(OUT_DIR, "read_tracking.tsv"))

message_time("DADA2 processing complete. Outputs written to ", normalizePath(OUT_DIR))
