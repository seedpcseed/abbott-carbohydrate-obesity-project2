# Simple phyloseq installation for radian
cat("=== Installing phyloseq for radian ===\n")

# Check current setup
cat("R version:", R.version.string, "\n")
cat("Library paths:", paste(.libPaths(), collapse=" | "), "\n")

# Install BiocManager if needed
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages('BiocManager', repos='https://cloud.r-project.org/')
}

# Try installing phyloseq directly without dependencies that are causing problems
cat("Attempting phyloseq installation...\n")
tryCatch({
  BiocManager::install('phyloseq', ask = FALSE, update = FALSE, dependencies = FALSE)
  cat("phyloseq installation completed.\n")
}, error = function(e) {
  cat("Direct installation failed, trying alternative approach...\n")
  cat("Error:", e$message, "\n")

  # Try installing just the R dependencies first
  cat("Installing core R packages...\n")
  install.packages(c('cluster', 'foreach', 'igraph', 'multtest', 'vegan'),
                   repos='https://cloud.r-project.org/', dependencies=FALSE)

  # Try Bioconductor packages one by one
  bioc_packages <- c('Biobase', 'BiocGenerics', 'biomformat', 'Biostrings', 'multtest')
  for(pkg in bioc_packages) {
    cat("Installing", pkg, "...\n")
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE, dependencies = FALSE)
    }, error = function(e2) {
      cat("Failed to install", pkg, ":", e2$message, "\n")
    })
  }

  # Finally try phyloseq again
  cat("Attempting phyloseq again...\n")
  BiocManager::install('phyloseq', ask = FALSE, update = FALSE, dependencies = FALSE)
})

# Check final status
cat("\n=== FINAL STATUS ===\n")
packages_to_check <- c("phyloseq", "ape", "vegan", "Biostrings")
for (pkg in packages_to_check) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste(pkg, "- INSTALLED\n"))
    } else {
        cat(paste(pkg, "- NOT FOUND\n"))
    }
}

cat("\nInstallation script completed.\n")