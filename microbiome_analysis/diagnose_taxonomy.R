#!/usr/bin/env Rscript

library(dada2)
library(Biostrings)

# Load the sequence table
cat("Loading sequence table...\n")
load("seqtab.nochim.RData")

# Check the structure of seqtab.nochim
cat("Sequence table dimensions:", dim(seqtab.nochim), "\n")
cat("Number of sequences:", ncol(seqtab.nochim), "\n")
cat("Number of samples:", nrow(seqtab.nochim), "\n")

# Check if sequences are empty or have issues
cat("Checking sequence lengths...\n")
seq_lengths <- nchar(colnames(seqtab.nochim))
cat("Sequence length summary:\n")
print(summary(seq_lengths))
cat("Any sequences with length 0:", any(seq_lengths == 0), "\n")
cat("Any sequences with length < 0:", any(seq_lengths < 0), "\n")

# Check for any NA or empty sequences
cat("Checking for NA or empty sequences...\n")
has_na <- any(is.na(colnames(seqtab.nochim)))
has_empty <- any(colnames(seqtab.nochim) == "")
cat("Has NA sequences:", has_na, "\n")
cat("Has empty sequences:", has_empty, "\n")

# Check the reference database
cat("Checking reference database...\n")
ref_file <- "references/silva_nr99_v138.2_toGenus_trainset.fa.gz"
if(file.exists(ref_file)) {
  cat("Reference file exists\n")
  cat("Reference file size:", file.size(ref_file), "bytes\n")
  
  # Try to read a few sequences from the reference
  tryCatch({
    ref_seqs <- readDNAStringSet(ref_file, nrec = 5)
    cat("Successfully read first 5 reference sequences\n")
    cat("Reference sequence lengths:", width(ref_seqs), "\n")
  }, error = function(e) {
    cat("Error reading reference sequences:", e$message, "\n")
  })
} else {
  cat("Reference file does not exist!\n")
}

# Try to identify problematic sequences
cat("Checking for problematic sequences...\n")
problematic_seqs <- which(seq_lengths <= 0 | is.na(seq_lengths))
if(length(problematic_seqs) > 0) {
  cat("Found", length(problematic_seqs), "problematic sequences at positions:", problematic_seqs, "\n")
  cat("These sequences:", colnames(seqtab.nochim)[problematic_seqs], "\n")
}

# Check if seqs variable exists and its structure
if(exists("seqs")) {
  cat("seqs variable exists with length:", length(seqs), "\n")
  cat("seqs class:", class(seqs), "\n")
} else {
  cat("seqs variable does not exist\n")
}

# Check if nonchimeric_seqs exists
if(exists("nonchimeric_seqs")) {
  cat("nonchimeric_seqs exists with length:", length(nonchimeric_seqs), "\n")
  cat("nonchimeric_seqs class:", class(nonchimeric_seqs), "\n")
} else {
  cat("nonchimeric_seqs variable does not exist\n")
}
