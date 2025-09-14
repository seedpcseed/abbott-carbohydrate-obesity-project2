# Final approach: Install phyloseq without problematic dependencies
cat("=== Final phyloseq installation attempt ===\n")

# Clean up any failed installations
unlink('/usr/lib/R/site-library/phyloseq', recursive=TRUE, force=TRUE)
unlink('/usr/lib/R/site-library/biomformat', recursive=TRUE, force=TRUE)
unlink('/usr/lib/R/site-library/rhdf5', recursive=TRUE, force=TRUE)

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos='https://cloud.r-project.org/')
}

# Strategy: Install phyloseq from an older version that doesn't require rhdf5/biomformat
# Or use alternative approach with pre-compiled binaries

# First, let's try installing from CRAN instead of Bioconductor
cat("Trying to install from alternative sources...\n")

# Install the minimal dependencies we know work
essential_deps <- c('ade4', 'cluster', 'foreach', 'igraph', 'plyr', 'reshape2', 'vegan')
cat("Installing essential dependencies...\n")
install.packages(essential_deps, repos='https://cloud.r-project.org/', dependencies=FALSE)

# Install Bioconductor essentials that don't require HDF5
bioc_deps <- c('Biobase', 'BiocGenerics', 'Biostrings', 'IRanges', 'S4Vectors')
cat("Installing Bioconductor dependencies...\n")
for(pkg in bioc_deps) {
    cat("Installing", pkg, "...\n")
    tryCatch({
        BiocManager::install(pkg, ask=FALSE, dependencies=FALSE, update=FALSE)
    }, error = function(e) {
        cat("Warning:", pkg, "failed -", e$message, "\n")
    })
}

# Try to manually create a minimal phyloseq-like environment
# This is a workaround if phyloseq itself can't be installed
cat("Checking if we have the core functionality needed...\n")

# Check what we have
required_for_microbiome <- c('vegan', 'ade4', 'Biostrings', 'cluster')
available <- sapply(required_for_microbiome, function(x) requireNamespace(x, quietly=TRUE))
cat("Available packages for microbiome analysis:\n")
print(available)

if(all(available)) {
    cat("\nGood news! You have the essential packages for microbiome analysis.\n")
    cat("You can proceed with most phyloseq-style analyses using these packages.\n")
} else {
    cat("Some essential packages are missing.\n")
}

# Final attempt at phyloseq
cat("\nFinal attempt at phyloseq installation...\n")
tryCatch({
    # Try with older bioconductor version or force install
    BiocManager::install('phyloseq', ask=FALSE, dependencies=FALSE, type='source')
}, error = function(e) {
    cat("Phyloseq installation failed. Using alternative approach.\n")

    # Create a simple replacement function
    cat("Creating basic phyloseq alternative functions...\n")

    # Write a small helper file
    writeLines('
# Basic phyloseq alternative functions
library(vegan)
library(Biostrings)

# Simple function to create phyloseq-like object
create_phyloseq_like <- function(otu_table, taxonomy, metadata=NULL) {
  result <- list(
    otu_table = otu_table,
    taxonomy = taxonomy,
    metadata = metadata
  )
  class(result) <- "phyloseq_like"
  return(result)
}

# Basic diversity calculations
calc_diversity <- function(otu_table, method="shannon") {
  diversity(otu_table, index=method)
}

cat("Basic phyloseq alternatives loaded. Use create_phyloseq_like() and calc_diversity()\\n")
', '/root/projects/abbott-carbohydrate-obesity-project2/microbiome_analysis/phyloseq_alternative.R')

    source('/root/projects/abbott-carbohydrate-obesity-project2/microbiome_analysis/phyloseq_alternative.R')
    cat("Phyloseq alternative functions created and sourced.\n")
})

cat("\n=== FINAL STATUS ===\n")
key_packages <- c("phyloseq", "vegan", "Biostrings", "ade4", "cluster")
for (pkg in key_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste(pkg, "- AVAILABLE\n"))
    } else {
        cat(paste(pkg, "- NOT AVAILABLE\n"))
    }
}