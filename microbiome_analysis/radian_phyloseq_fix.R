# Radian phyloseq installation fix
# Run this in your radian session

# Check if essential packages are available
essential_packages <- c("vegan", "Biostrings", "ade4", "cluster", "reshape2")
missing_packages <- essential_packages[!sapply(essential_packages, requireNamespace, quietly=TRUE)]

if(length(missing_packages) > 0) {
    cat("Installing missing essential packages:", paste(missing_packages, collapse=", "), "\n")
    install.packages(missing_packages, repos="https://cloud.r-project.org/")
}

# Try phyloseq installation without dependencies that cause problems
if(!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org/")
}

# Set up environment to avoid HDF5 issues
Sys.setenv("HDF5_USE_FILE_LOCKING"="FALSE")

# Try installing phyloseq with minimal dependencies
cat("Attempting phyloseq installation without problematic dependencies...\n")
BiocManager::install("phyloseq", ask=FALSE, dependencies=FALSE, update=FALSE)

# Check if it worked
if(requireNamespace("phyloseq", quietly=TRUE)) {
    cat("SUCCESS: phyloseq is now available!\n")
    cat("You can now use: library(phyloseq)\n")
} else {
    cat("phyloseq installation failed. Setting up alternative microbiome analysis environment...\n")

    # Load essential packages for microbiome analysis
    library(vegan)
    library(Biostrings)
    library(ade4)
    library(cluster)

    cat("Alternative microbiome analysis packages loaded:\n")
    cat("- vegan: for diversity analysis and ordination\n")
    cat("- Biostrings: for sequence handling\n")
    cat("- ade4: for multivariate analysis\n")
    cat("- cluster: for clustering analysis\n")

    cat("\nYou can perform most microbiome analyses with these packages.\n")
    cat("For phyloseq-style workflows, you can use vegan and Biostrings together.\n")
}

cat("\nSetup complete!\n")