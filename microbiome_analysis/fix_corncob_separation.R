#!/usr/bin/env Rscript

# Fix corncob separation issues
library(phyloseq)
library(corncob)
library(dplyr)

# Load the phyloseq object
cat("Loading phyloseq object...\n")
ps_raw <- readRDS("outputs/ps_raw.rds")

# Add SubjectID if missing
sdat <- sample_data(ps_raw) %>% as.data.frame()
if (!"SubjectID" %in% colnames(sdat) && "subjectID" %in% colnames(sdat)) {
  sdat$SubjectID <- as.factor(sdat$subjectID)
  sample_data(ps_raw) <- sample_data(sdat)
  cat("Created SubjectID from subjectID\n")
}

# Work at Genus level
ps_genus_counts <- tax_glom(ps_raw, taxrank = "Genus", NArm = FALSE)

# Apply more aggressive prevalence filtering to avoid separation
prev_threshold <- 0.15  # Increased from 0.05
min_prev <- ceiling(prev_threshold * nsamples(ps_genus_counts))
keep_taxa <- taxa_sums(ps_genus_counts) > 0 & rowSums(otu_table(ps_genus_counts) > 0) >= min_prev
ps_genus_filt <- prune_taxa(keep_taxa, ps_genus_counts)

cat("Original genus taxa:", ntaxa(ps_genus_counts), "\n")
cat("Filtered genus taxa:", ntaxa(ps_genus_filt), "\n")

# Function to check for separation issues
check_separation <- function(ps, formula_str) {
  sdat <- sample_data(ps) %>% as.data.frame()
  otu_mat <- otu_table(ps) %>% as.matrix()
  
  # Parse formula to get variables
  vars <- gsub("~", "", formula_str)
  vars <- trimws(unlist(strsplit(vars, "\\+")))
  vars <- trimws(gsub("\\(.*\\)", "", vars))  # Remove random effects
  vars <- unique(vars[vars != ""])
  
  problematic_taxa <- c()
  
  for(i in 1:ntaxa(ps)) {
    taxon_name <- taxa_names(ps)[i]
    abundances <- otu_mat[i, , drop = TRUE]
    
    # Check for separation in each variable
    for(var in vars) {
      if(var %in% colnames(sdat)) {
        var_levels <- levels(as.factor(sdat[[var]]))
        
        for(level in var_levels) {
          level_samples <- which(sdat[[var]] == level)
          if(length(level_samples) > 0 && all(level_samples <= length(abundances))) {
            level_abundances <- abundances[level_samples]
            
            # Check if all abundances are zero in this level
            if(all(level_abundances == 0)) {
              problematic_taxa <- c(problematic_taxa, taxon_name)
              cat("Separation detected: Taxon", taxon_name, "is all zeros in", var, "=", level, "\n")
              break
            }
          }
        }
      }
    }
  }
  
  return(unique(problematic_taxa))
}

# Check for separation issues
formula_str <- "~ group + timepoint + carbohydrate + SubjectID"
cat("\nChecking for separation issues...\n")
problematic_taxa <- check_separation(ps_genus_filt, formula_str)

if(length(problematic_taxa) > 0) {
  cat("Found", length(problematic_taxa), "problematic taxa with separation issues\n")
  cat("Removing problematic taxa...\n")
  ps_genus_clean <- prune_taxa(!taxa_names(ps_genus_filt) %in% problematic_taxa, ps_genus_filt)
  cat("Cleaned dataset has", ntaxa(ps_genus_clean), "taxa\n")
} else {
  cat("No separation issues detected\n")
  ps_genus_clean <- ps_genus_filt
}

# Test corncob with cleaned data
cat("\nTesting corncob with cleaned data...\n")

# Function to run corncob with error handling
run_corncob_safe <- function(ps, effect_var, fdr_cutoff = 0.1) {
  cat("Testing effect:", effect_var, "\n")
  
  # Create formulas
  all_vars <- c("group", "timepoint", "carbohydrate", "SubjectID")
  full_vars <- all_vars
  red_vars <- setdiff(all_vars, effect_var)
  
  full <- as.formula(paste("~", paste(full_vars, collapse = " + ")))
  red <- as.formula(paste("~", paste(red_vars, collapse = " + ")))
  
  # Try with different settings
  settings <- list(
    list(boot = FALSE, B = 0, fdr_cutoff = fdr_cutoff),
    list(boot = TRUE, B = 100, fdr_cutoff = fdr_cutoff),
    list(boot = FALSE, B = 0, fdr_cutoff = fdr_cutoff, method = "LRT")
  )
  
  for(i in seq_along(settings)) {
    cat("  Trying setting", i, "...\n")
    setting <- settings[[i]]
    
    tryCatch({
      res <- differentialTest(
        formula = full,
        phi.formula = full,
        formula_null = red,
        phi.formula_null = red,
        test = "Wald",
        data = ps,
        boot = setting$boot,
        B = setting$B,
        fdr_cutoff = setting$fdr_cutoff
      )
      
      cat("  SUCCESS with setting", i, "\n")
      return(tibble(
        taxon = res$significant_taxa,
        effect = effect_var,
        method = paste0("setting_", i)
      ))
    }, error = function(e) {
      cat("  FAILED with setting", i, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  cat("  All settings failed for", effect_var, "\n")
  return(tibble(taxon = character(0), effect = effect_var, method = "failed"))
}

# Test each effect
effects <- c("group", "timepoint", "carbohydrate")
results <- list()

for(effect in effects) {
  results[[effect]] <- run_corncob_safe(ps_genus_clean, effect)
}

# Combine results
all_results <- bind_rows(results)

if(nrow(all_results) > 0) {
  cat("\nCorncob analysis completed successfully!\n")
  cat("Found", nrow(all_results), "significant taxa across all effects\n")
  
  # Save results
  write_csv(all_results, "outputs/corncob_genus_hits_fixed.csv")
  cat("Results saved to outputs/corncob_genus_hits_fixed.csv\n")
  
  # Show summary
  cat("\nSummary by effect:\n")
  print(table(all_results$effect))
} else {
  cat("\nNo significant taxa found with any method\n")
}

# Alternative: Try with even more aggressive filtering
if(nrow(all_results) == 0) {
  cat("\nTrying with more aggressive filtering...\n")
  
  # Increase prevalence threshold
  prev_threshold <- 0.25
  min_prev <- ceiling(prev_threshold * nsamples(ps_genus_counts))
  keep_taxa <- taxa_sums(ps_genus_counts) > 0 & rowSums(otu_table(ps_genus_counts) > 0) >= min_prev
  ps_genus_aggressive <- prune_taxa(keep_taxa, ps_genus_counts)
  
  cat("Aggressively filtered dataset has", ntaxa(ps_genus_aggressive), "taxa\n")
  
  # Test with simplified model
  tryCatch({
    res_simple <- differentialTest(
      formula = ~ group,
      phi.formula = ~ group,
      formula_null = ~ 1,
      phi.formula_null = ~ 1,
      test = "Wald",
      data = ps_genus_aggressive,
      boot = FALSE,
      B = 0
    )
    
    cat("Simple model (group only) succeeded with", length(res_simple$significant_taxa), "significant taxa\n")
    
    simple_results <- tibble(
      taxon = res_simple$significant_taxa,
      effect = "group",
      method = "simple"
    )
    
    write_csv(simple_results, "outputs/corncob_genus_hits_simple.csv")
    cat("Simple results saved to outputs/corncob_genus_hits_simple.csv\n")
    
  }, error = function(e) {
    cat("Even simple model failed:", e$message, "\n")
  })
}

cat("\nDebugging complete. Check the output files for results.\n")
