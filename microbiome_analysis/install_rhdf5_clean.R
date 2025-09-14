# Alternative approach: Install without rhdf5 dependency
cat("Installing BiocManager if needed...\n")
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos='https://cloud.r-project.org/')
}

# Try to install phyloseq without rhdf5 first
cat("Installing ape...\n")
install.packages('ape', repos='https://cloud.r-project.org/')

# Install other phyloseq dependencies without rhdf5
cat("Installing phyloseq without rhdf5 dependency...\n")
BiocManager::install(c('Biostrings', 'biomformat'), ask=FALSE)

# Now try phyloseq
cat("Installing phyloseq...\n")
BiocManager::install('phyloseq', ask=FALSE)

# Verify what we have
packages_to_check <- c("ape", "Biostrings", "biomformat", "phyloseq")
cat("\n=== INSTALLATION STATUS ===\n")
for (pkg in packages_to_check) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste(pkg, "- INSTALLED\n"))
    } else {
        cat(paste(pkg, "- FAILED\n"))
    }
}

cat("\n=== Alternative: Try rhdf5 separately ===\n")
# Simple rhdf5 installation
BiocManager::install('rhdf5', ask=FALSE)
if (requireNamespace("rhdf5", quietly = TRUE)) {
    cat("rhdf5 - INSTALLED\n")
} else {
    cat("rhdf5 - FAILED (this is optional for most phyloseq functions)\n")
}

