#!/usr/bin/env Rscript

# Simple fix for corncob separation issues
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

# Apply aggressive prevalence filtering to avoid separation
prev_threshold <- 0.2  # Very high threshold
min_prev <- ceiling(prev_threshold * nsamples(ps_genus_counts))
keep_taxa <- taxa_sums(ps_genus_counts) > 0 & rowSums(otu_table(ps_genus_counts) > 0) >= min_prev
ps_genus_filt <- prune_taxa(keep_taxa, ps_genus_counts)

cat("Original genus taxa:", ntaxa(ps_genus_counts), "\n")
cat("Aggressively filtered genus taxa:", ntaxa(ps_genus_filt), "\n")

# Test with a very simple model first
cat("\nTesting with simple group-only model...\n")
tryCatch({
  res_simple <- differentialTest(
    formula = ~ group,
    phi.formula = ~ group,
    formula_null = ~ 1,
    phi.formula_null = ~ 1,
    test = "Wald",
    data = ps_genus_filt,
    boot = FALSE,
    B = 0,
    fdr_cutoff = 0.1
  )
  
  cat("Simple model SUCCESS!\n")
  cat("Found", length(res_simple$significant_taxa), "significant taxa\n")
  
  if(length(res_simple$significant_taxa) > 0) {
    simple_results <- tibble(
      taxon = res_simple$significant_taxa,
      effect = "group",
      method = "simple_group_only"
    )
    
    write_csv(simple_results, "outputs/corncob_genus_hits_simple.csv")
    cat("Results saved to outputs/corncob_genus_hits_simple.csv\n")
    print(simple_results)
  }
  
}, error = function(e) {
  cat("Simple model FAILED:", e$message, "\n")
})

# Test with timepoint
cat("\nTesting with timepoint model...\n")
tryCatch({
  res_time <- differentialTest(
    formula = ~ timepoint,
    phi.formula = ~ timepoint,
    formula_null = ~ 1,
    phi.formula_null = ~ 1,
    test = "Wald",
    data = ps_genus_filt,
    boot = FALSE,
    B = 0,
    fdr_cutoff = 0.1
  )
  
  cat("Timepoint model SUCCESS!\n")
  cat("Found", length(res_time$significant_taxa), "significant taxa\n")
  
  if(length(res_time$significant_taxa) > 0) {
    time_results <- tibble(
      taxon = res_time$significant_taxa,
      effect = "timepoint",
      method = "simple_timepoint_only"
    )
    
    write_csv(time_results, "outputs/corncob_genus_hits_timepoint.csv")
    cat("Results saved to outputs/corncob_genus_hits_timepoint.csv\n")
    print(time_results)
  }
  
}, error = function(e) {
  cat("Timepoint model FAILED:", e$message, "\n")
})

# Test with carbohydrate
cat("\nTesting with carbohydrate model...\n")
tryCatch({
  res_carb <- differentialTest(
    formula = ~ carbohydrate,
    phi.formula = ~ carbohydrate,
    formula_null = ~ 1,
    phi.formula_null = ~ 1,
    test = "Wald",
    data = ps_genus_filt,
    boot = FALSE,
    B = 0,
    fdr_cutoff = 0.1
  )
  
  cat("Carbohydrate model SUCCESS!\n")
  cat("Found", length(res_carb$significant_taxa), "significant taxa\n")
  
  if(length(res_carb$significant_taxa) > 0) {
    carb_results <- tibble(
      taxon = res_carb$significant_taxa,
      effect = "carbohydrate",
      method = "simple_carbohydrate_only"
    )
    
    write_csv(carb_results, "outputs/corncob_genus_hits_carbohydrate.csv")
    cat("Results saved to outputs/corncob_genus_hits_carbohydrate.csv\n")
    print(carb_results)
  }
  
}, error = function(e) {
  cat("Carbohydrate model FAILED:", e$message, "\n")
})

# Now try with multiple variables but simpler structure
cat("\nTesting with group + timepoint model...\n")
tryCatch({
  res_gt <- differentialTest(
    formula = ~ group + timepoint,
    phi.formula = ~ group + timepoint,
    formula_null = ~ group,
    phi.formula_null = ~ group,
    test = "Wald",
    data = ps_genus_filt,
    boot = FALSE,
    B = 0,
    fdr_cutoff = 0.1
  )
  
  cat("Group + timepoint model SUCCESS!\n")
  cat("Found", length(res_gt$significant_taxa), "significant taxa\n")
  
  if(length(res_gt$significant_taxa) > 0) {
    gt_results <- tibble(
      taxon = res_gt$significant_taxa,
      effect = "timepoint_adjusted_group",
      method = "group_timepoint"
    )
    
    write_csv(gt_results, "outputs/corncob_genus_hits_group_timepoint.csv")
    cat("Results saved to outputs/corncob_genus_hits_group_timepoint.csv\n")
    print(gt_results)
  }
  
}, error = function(e) {
  cat("Group + timepoint model FAILED:", e$message, "\n")
})

cat("\nSimple corncob testing complete!\n")
cat("Check the output files for results from each model.\n")
