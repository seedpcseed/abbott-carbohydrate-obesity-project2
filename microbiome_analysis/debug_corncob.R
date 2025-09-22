#!/usr/bin/env Rscript

# Debug script for corncob differentialTest error
library(phyloseq)
library(corncob)
library(dplyr)

# Load the phyloseq object
cat("Loading phyloseq object...\n")
ps_raw <- readRDS("outputs/ps_raw.rds")

# Check the data structure
cat("=== PHYLOSEQ OBJECT STRUCTURE ===\n")
cat("Number of samples:", nsamples(ps_raw), "\n")
cat("Number of taxa:", ntaxa(ps_raw), "\n")
cat("Sample variables:", paste(sample_variables(ps_raw), collapse = ", "), "\n")

# Check for missing values
sdat <- sample_data(ps_raw) %>% as.data.frame()
cat("\n=== MISSING VALUES ===\n")
for(var in c("group", "timepoint", "carbohydrate", "SubjectID")) {
  if(var %in% colnames(sdat)) {
    missing_count <- sum(is.na(sdat[[var]]))
    cat(var, ":", missing_count, "missing values\n")
  } else {
    cat(var, ": NOT FOUND\n")
  }
}

# Check factor levels
cat("\n=== FACTOR LEVELS ===\n")
for(var in c("group", "timepoint", "carbohydrate", "SubjectID")) {
  if(var %in% colnames(sdat)) {
    cat(var, "levels:", paste(levels(as.factor(sdat[[var]])), collapse = ", "), "\n")
  }
}

# Create SubjectID if missing
if (!"SubjectID" %in% colnames(sdat) && "subjectID" %in% colnames(sdat)) {
  sdat$SubjectID <- as.factor(sdat$subjectID)
  sample_data(ps_raw) <- sample_data(sdat)
  cat("Created SubjectID from subjectID\n")
}

# Work at Genus level
cat("\n=== GENUS AGGREGATION ===\n")
ps_genus_counts <- tax_glom(ps_raw, taxrank = "Genus", NArm = TRUE)
cat("Genus taxa count:", ntaxa(ps_genus_counts), "\n")

# Apply prevalence filtering
prev_threshold <- 0.05
min_prev <- ceiling(prev_threshold * nsamples(ps_genus_counts))
keep_taxa <- taxa_sums(ps_genus_counts) > 0 & rowSums(otu_table(ps_genus_counts) > 0) >= min_prev
ps_genus_filt <- prune_taxa(keep_taxa, ps_genus_counts)
cat("Filtered genus taxa count:", ntaxa(ps_genus_filt), "\n")

# Check for zero-sum taxa
otu_mat <- otu_table(ps_genus_filt) %>% as.matrix()
zero_sum_taxa <- rowSums(otu_mat) == 0
cat("Taxa with zero total abundance:", sum(zero_sum_taxa), "\n")

# Check for taxa with very low variance
taxa_var <- apply(otu_mat, 1, function(x) var(x))
low_var_taxa <- taxa_var < 1e-10
cat("Taxa with very low variance:", sum(low_var_taxa), "\n")

# Test with a single taxon first
cat("\n=== TESTING WITH SINGLE TAXON ===\n")
single_taxon <- names(which(rowSums(otu_mat) > 0))[1]
cat("Testing with taxon:", single_taxon, "\n")

# Subset to single taxon
ps_single <- prune_taxa(single_taxon, ps_genus_filt)

# Try a simple model first
tryCatch({
  cat("Testing simple model...\n")
  res_simple <- differentialTest(
    formula = ~ group,
    phi.formula = ~ group,
    formula_null = ~ 1,
    phi.formula_null = ~ 1,
    test = "Wald",
    data = ps_single,
    boot = FALSE,
    B = 0
  )
  cat("Simple model SUCCESS\n")
}, error = function(e) {
  cat("Simple model FAILED:", e$message, "\n")
})

# Try with more complex model
ps_genus_filt_sdcrdc <- subset_samples(ps_genus_filt, carbohydrate%in%c("SDC", "RDC"))
tryCatch({
  cat("Testing complex model...\n")
  res_complex <- differentialTest(
    formula = ~ group + timepoint + carbohydrate,
    phi.formula = ~ group + timepoint + carbohydrate,
    formula_null = ~ timepoint, 
    phi.formula_null = ~ group + timepoint + carbohydrate,
    test = "Wald",
    data = ps_genus_filt_sdcrdc,
    boot = FALSE,
    B = 0
  )
  cat("Complex model SUCCESS\n")
}, error = function(e) {
  cat("Complex model FAILED:", e$message, "\n")
})

# Check for problematic taxa patterns
cat("\n=== CHECKING FOR PROBLEMATIC PATTERNS ===\n")

# Check for taxa that are all zeros in one group
sdat_filt <- sample_data(ps_genus_filt) %>% as.data.frame()
group_levels <- levels(as.factor(sdat_filt$group))

problematic_taxa <- c()
for(i in 1:ntaxa(ps_genus_filt)) {
  taxon_name <- taxa_names(ps_genus_filt)[i]
  abundances <- otu_mat[i, ]
  
  # Check if taxon is all zeros in any group
  for(grp in group_levels) {
    grp_samples <- which(sdat_filt$group == grp)
    if(length(grp_samples) > 0) {
      grp_abundances <- abundances[grp_samples]
      
      if(all(grp_abundances == 0)) {
        problematic_taxa <- c(problematic_taxa, taxon_name)
        cat("Taxon", taxon_name, "is all zeros in group", grp, "\n")
      }
    }
  }
}

cat("Total problematic taxa found:", length(unique(problematic_taxa)), "\n")

# Suggest solution
cat("\n=== SUGGESTED SOLUTION ===\n")
if(length(unique(problematic_taxa)) > 0) {
  cat("Remove problematic taxa and retry:\n")
  cat("ps_genus_filt_clean <- prune_taxa(!taxa_names(ps_genus_filt) %in% unique(problematic_taxa), ps_genus_filt)\n")
} else {
  cat("No obviously problematic taxa found. Try:\n")
  cat("1. Check for convergence issues with boot = TRUE\n")
  cat("2. Use fewer taxa (increase prevalence threshold)\n")
  cat("3. Simplify the model\n")
}
