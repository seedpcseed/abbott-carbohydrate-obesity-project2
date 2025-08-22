#!/usr/bin/env Rscript

library(dada2)
library(Biostrings)

# Load the original sequence table
cat("Loading original sequence table...\n")
load("seqtab.RData")

# Read the non-chimeric sequences
cat("Reading non-chimeric sequences...\n")
nonchimeric_seqs <- readDNAStringSet("repseqs.nochimera.fasta")

# Find which columns in seqtab correspond to non-chimeric sequences
cat("Finding non-chimeric sequences in sequence table...\n")
w <- which(colnames(seqtab) %in% as.character(nonchimeric_seqs))
cat("Found", length(w), "non-chimeric sequences out of", ncol(seqtab), "total sequences\n")

# Create the non-chimeric sequence table
seqtab.nochim <- seqtab[, w]

# Get the sequences
seqs <- colnames(seqtab.nochim)

# Save the regenerated data
cat("Saving regenerated sequence table...\n")
save(mergers, seqtab, seqtab.nochim, seqs, nonchimeric_seqs, file = "seqtab.nochim.RData")

# Check the structure
cat("Sequence table dimensions:", dim(seqtab.nochim), "\n")
cat("Number of sequences:", ncol(seqtab.nochim), "\n")
cat("Number of samples:", nrow(seqtab.nochim), "\n")

# Check sequence lengths
seq_lengths <- nchar(colnames(seqtab.nochim))
cat("Sequence length summary:\n")
print(summary(seq_lengths))
cat("Any sequences with length 0:", any(seq_lengths == 0), "\n")
cat("Any sequences with length < 0:", any(seq_lengths < 0), "\n")

cat("Regeneration complete!\n")
