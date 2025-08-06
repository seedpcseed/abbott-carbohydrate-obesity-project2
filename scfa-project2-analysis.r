library(tidyverse)
library(readxl)
library(janitor)
library(data.table)
library(ggpubr)
library(rstatix)
library(viridis)
library(RColorBrewer)
library(lme4)
library(lmerTest)

# ===============================
# Read in metadata
# ===============================

map <- read_excel("metadata/SEED LAB HMM Submission Sheet 202405.xlsx", skip = 9) %>% clean_names()
# replace column names with cleaned names
colnames(map) <- c("sampleid", "group", "bile_acids_panel", "pfbbr_scfa_panel", "tryp_indole_panel", "tms_panel_1", "tms_panel_2", "tms_panel_3", "untargeted_panel", "custom_targeted_panel", "low_abundant_metabolite_panel", "sample_origin", "sample_type", "mouse_number", "services_requested_per_sample", "mouse_strain", "mouse_flora", "volume_ul", "sample_weight_mg", "experiment_treatment")
# remove rows where experiment_sample_name is NA
map <- map %>% filter(!is.na(sampleid))
# remove rows where experiment_sample_name is empty
map <- map %>% filter(sampleid != "")
#remove core lab information
map <- map[-(1:3), ]

#split out sample name information



# Extract information from sample names
map <- map %>%
  mutate(
    # Extract subject ID (digits and letter before first space)
    subject = str_extract(sampleid, "^[0-9]+[a-zA-Z]*"),
    
    # Extract the part between spaces (e.g., "r1", "s2", "nc")
    carb_repeat = str_extract(sampleid, "(?<= )[rs]\\d|nc(?= )"),
    
    # Create carbohydrate_type column
    carbohydrate_type = case_when(
      str_detect(carb_repeat, "^r") ~ "rapid_digestible",
      str_detect(carb_repeat, "^s") ~ "slow_digestible", 
      str_detect(carb_repeat, "^nc") ~ "no_carbohydrate",
      TRUE ~ NA_character_
    ),
    
    # Extract technical repeat number (handle "nc" case where there's no number)
    technical_repeat = case_when(
      str_detect(carb_repeat, "^nc") ~ 1L,  # nc gets repeat 1
      TRUE ~ as.integer(str_extract(carb_repeat, "\\d+"))
    ),
    
    # Extract timepoint (after second space, before "h")
    timepoint_hr = as.integer(str_extract(sampleid, "(?<= )\\d+(?=h$)"))
  ) %>%
  # Remove the temporary carb_repeat column
  select(-carb_repeat)

map <- map %>% select(sampleid, group, subject, carbohydrate_type, technical_repeat, timepoint_hr)

map2 <- read_excel("metadata/PS2720_xtra_samples_HMM Submission Sheet.xlsx", skip = 9) %>% clean_names()
# replace column names with cleaned names
colnames(map2) <- c("sampleid", "group", "bile_acids_panel", "pfbbr_scfa_panel", "tryp_indole_panel", "tms_panel_1", "tms_panel_2", "tms_panel_3", "untargeted_panel", "custom_targeted_panel", "low_abundant_metabolite_panel", "sample_origin", "sample_type", "mouse_number", "services_requested_per_sample", "mouse_strain", "mouse_flora", "volume_ul", "sample_weight_mg", "experiment_treatment")
# remove rows where experiment_sample_name is NA
map2 <- map2 %>% filter(!is.na(sampleid))
# remove rows where experiment_sample_name is empty
map2 <- map2 %>% filter(sampleid != "")
#remove core lab information
map2 <- map2[-(1:3), ]

#split out sample name information
map2 <- map2 %>%
  mutate(
    # Extract subject ID (digits and letter before first space)
    subject = str_extract(sampleid, "^[0-9]+[a-zA-Z]*"),
    
    # Extract the part between spaces (e.g., "r1", "s2", "nc")
    carb_repeat = str_extract(sampleid, "(?<= )[rs]\\d|nc(?= )"),
    
    # Create carbohydrate_type column
    carbohydrate_type = case_when(
      str_detect(carb_repeat, "^r") ~ "rapid_digestible",
      str_detect(carb_repeat, "^s") ~ "slow_digestible", 
      str_detect(carb_repeat, "^nc") ~ "no_carbohydrate",
      TRUE ~ NA_character_
    ),
    
    # Extract technical repeat number (handle "nc" case where there's no number)
    technical_repeat = case_when(
      str_detect(carb_repeat, "^nc") ~ 1L,  # nc gets repeat 1
      TRUE ~ as.integer(str_extract(carb_repeat, "\\d+"))
    ),
    
    # Extract timepoint (after second space, before "h")
    timepoint_hr = as.integer(str_extract(sampleid, "(?<= )\\d+(?=h$)"))
  ) %>%
  # Remove the temporary carb_repeat column
  select(-carb_repeat)

map2 <- map2 %>% select(sampleid, group, subject, carbohydrate_type, technical_repeat, timepoint_hr)

map <- bind_rows(map, map2)

# Make group names lowercase for consistency
map <- map %>%
  mutate(group = tolower(group))

# ===============================
# Read in data
# ===============================
data <- fread("data/removed_qcs_quant_results_20250722_PFBBr_PS2720_20250731.csv")

sample_data <- left_join(map, data, by = "sampleid")

# Ensure group names are lowercase
sample_data <- sample_data %>%
  mutate(group = tolower(group))


# ===============================
# Calculate summary statistics
# ===============================

# Define SCFA analytes (excluding non-SCFA compounds)
scfa_analytes <- c("acetate", "butyrate", "propionate", "5aminovalerate", "succinate")

# Convert data to long format and average technical replicates
sample_data_long <- sample_data %>%
  pivot_longer(cols = all_of(scfa_analytes), 
               names_to = "analyte", 
               values_to = "concentration") %>%
  filter(!is.na(concentration)) %>%
  # Average technical replicates FIRST before any statistical analysis
  group_by(subject, group, carbohydrate_type, timepoint_hr, analyte) %>%
  summarise(concentration = mean(concentration, na.rm = TRUE), .groups = "drop")

# Summary statistics by group
summary_by_group <- sample_data_long %>%
  group_by(group, analyte) %>%
  summarise(
    n = n(),
    mean = mean(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    sem = sd / sqrt(n),
    q25 = quantile(concentration, 0.25, na.rm = TRUE),
    q75 = quantile(concentration, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Order so control and case are adjacent for each analyte
  arrange(analyte, group)

# Summary statistics by carbohydrate type
summary_by_carb <- sample_data_long %>%
  group_by(carbohydrate_type, analyte) %>%
  summarise(
    n = n(),
    mean = mean(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    sem = sd / sqrt(n),
    q25 = quantile(concentration, 0.25, na.rm = TRUE),
    q75 = quantile(concentration, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Order so the 3 carbohydrate types are in series with each analyte
  arrange(analyte, carbohydrate_type)

# Summary statistics by timepoint
summary_by_time <- sample_data_long %>%
  group_by(timepoint_hr, analyte) %>%
  summarise(
    n = n(),
    mean = mean(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    sem = sd / sqrt(n),
    q25 = quantile(concentration, 0.25, na.rm = TRUE),
    q75 = quantile(concentration, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Order so that 0h and 48h are adjacent per analyte
  arrange(analyte, timepoint_hr)

# Combined summary statistics by all three factors
summary_combined <- sample_data_long %>%
  group_by(group, carbohydrate_type, timepoint_hr, analyte) %>%
  summarise(
    n = n(),
    mean = mean(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    sd = sd(concentration, na.rm = TRUE),
    sem = sd / sqrt(n),
    q25 = quantile(concentration, 0.25, na.rm = TRUE),
    q75 = quantile(concentration, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# ===============================
# Statistical Analysis
# ===============================

# Prepare data with proper factor levels for statistical testing
sample_data_long <- sample_data_long %>%
  mutate(
    group = factor(group, levels = c("control", "case")),
    carbohydrate_type = factor(carbohydrate_type, 
                              levels = c("no_carbohydrate", "rapid_digestible", "slow_digestible")),
    timepoint_hr = factor(timepoint_hr),
    subject = factor(subject)  # Add subject as factor for random effects
  )

# Statistical tests by group (t-tests)
group_stats <- sample_data_long %>%
  group_by(analyte) %>%
  t_test(concentration ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Statistical tests by carbohydrate type (ANOVA with post-hoc)
carb_stats <- sample_data_long %>%
  group_by(analyte) %>%
  anova_test(concentration ~ carbohydrate_type) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Post-hoc comparisons for carbohydrate type (vs no_carbohydrate)
carb_posthoc <- sample_data_long %>%
  group_by(analyte) %>%
  pairwise_t_test(concentration ~ carbohydrate_type, 
                  ref.group = "no_carbohydrate") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Three-way ANOVA: Group x Carbohydrate x Time interaction
interaction_stats <- sample_data_long %>%
  group_by(analyte) %>%
  anova_test(concentration ~ group * carbohydrate_type * timepoint_hr) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Mixed-effects models accounting for subject as random effect
mixed_effects_results <- list()

for (analyte_name in scfa_analytes) {
  analyte_data <- sample_data_long %>% filter(analyte == analyte_name)
  
  # Full model with all interactions and subject as random intercept
  full_model <- lmer(concentration ~ group * carbohydrate_type * timepoint_hr + 
                     (1 | subject), 
                     data = analyte_data, 
                     REML = FALSE)
  
  # Model without three-way interaction
  no_three_way <- lmer(concentration ~ group + carbohydrate_type + timepoint_hr +
                       group:carbohydrate_type + group:timepoint_hr + 
                       carbohydrate_type:timepoint_hr + (1 | subject), 
                       data = analyte_data, 
                       REML = FALSE)
  
  # Model with only main effects
  main_effects <- lmer(concentration ~ group + carbohydrate_type + timepoint_hr + 
                       (1 | subject), 
                       data = analyte_data, 
                       REML = FALSE)
  
  # Store results
  mixed_effects_results[[analyte_name]] <- list(
    full_model = full_model,
    no_three_way = no_three_way,
    main_effects = main_effects,
    anova_full = anova(full_model),
    summary_full = summary(full_model)
  )
}

# Subject-level summary statistics
subject_summary <- sample_data_long %>%
  group_by(subject, group, analyte) %>%
  summarise(
    mean_conc = mean(concentration, na.rm = TRUE),
    n_observations = n(),
    .groups = "drop"
  ) %>%
  group_by(group, analyte) %>%
  summarise(
    n_subjects = n(),
    mean_subject_means = mean(mean_conc, na.rm = TRUE),
    sd_subject_means = sd(mean_conc, na.rm = TRUE),
    sem_subjects = sd_subject_means / sqrt(n_subjects),
    .groups = "drop"
  )

# Within-subject changes (baseline to 48h)
within_subject_changes <- sample_data_long %>%
  filter(timepoint_hr %in% c("0", "48")) %>%
  mutate(timepoint_hr = as.character(timepoint_hr)) %>%  # Convert back to character for pivot
  select(subject, group, carbohydrate_type, timepoint_hr, analyte, concentration) %>%
  pivot_wider(names_from = timepoint_hr, 
              values_from = concentration, 
              names_prefix = "time_") %>%
  mutate(change_0_to_48 = time_48 - time_0) %>%
  filter(!is.na(change_0_to_48))

# Statistical tests for within-subject changes
within_subject_stats <- within_subject_changes %>%
  group_by(group, carbohydrate_type, analyte) %>%
  summarise(
    n = n(),
    mean_change = mean(change_0_to_48, na.rm = TRUE),
    sd_change = sd(change_0_to_48, na.rm = TRUE),
    sem_change = sd_change / sqrt(n),
    t_stat = mean_change / sem_change,
    p_value = 2 * pt(-abs(t_stat), df = n - 1),  # two-tailed t-test
    .groups = "drop"
  ) %>%
  group_by(analyte) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(significance = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**", 
    p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Create results directory and save combined summary
dir.create("results", showWarnings = FALSE)
write.csv(summary_combined, "results/combined_summary_statistics.csv", row.names = FALSE)

# Print statistical results
cat("\n=== STATISTICAL ANALYSIS RESULTS ===\n")
cat("\nGroup Comparisons (Control vs Case):\n")
print(group_stats)
cat("\nCarbohydrate Type ANOVA:\n")
print(carb_stats)
cat("\nCarbohydrate Post-hoc (vs No Carbohydrate):\n")
print(carb_posthoc)
cat("\nThree-way Interaction Analysis:\n")
interaction_display <- as.data.frame(interaction_stats) %>%
  select(analyte, Effect, p, p.adj, p.adj.signif)
print(interaction_display)

# Case-only temporal analysis (0h vs 48h)
case_only_changes <- sample_data_long %>%
  filter(group == "case", timepoint_hr %in% c("0", "48")) %>%
  mutate(timepoint_hr = as.character(timepoint_hr)) %>%
  pivot_wider(names_from = timepoint_hr, 
              values_from = concentration, 
              names_prefix = "time_") %>%
  mutate(change_0_to_48 = time_48 - time_0) %>%
  filter(!is.na(change_0_to_48))

# Statistical tests for case-only temporal changes
case_only_stats <- case_only_changes %>%
  group_by(carbohydrate_type, analyte) %>%
  summarise(
    n = n(),
    mean_change = mean(change_0_to_48, na.rm = TRUE),
    sd_change = sd(change_0_to_48, na.rm = TRUE),
    sem_change = sd_change / sqrt(n),
    t_stat = mean_change / sem_change,
    p_value = 2 * pt(-abs(t_stat), df = n - 1),  # two-tailed t-test
    .groups = "drop"
  ) %>%
  group_by(analyte) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(significance = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**", 
    p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Paired t-tests for case-only changes by analyte (pooled across carb types)
case_only_pooled_stats <- case_only_changes %>%
  group_by(analyte) %>%
  summarise(
    n = n(),
    mean_change = mean(change_0_to_48, na.rm = TRUE),
    sd_change = sd(change_0_to_48, na.rm = TRUE),
    sem_change = sd_change / sqrt(n),
    t_stat = mean_change / sem_change,
    p_value = 2 * pt(-abs(t_stat), df = n - 1),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(significance = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**", 
    p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Mixed-effects models for case-only temporal analysis
case_only_mixed_effects <- list()

# Filter data to case group only
case_only_data <- sample_data_long %>%
  filter(group == "case") %>%
  mutate(
    timepoint_hr = factor(timepoint_hr),
    subject = factor(subject),
    carbohydrate_type = factor(carbohydrate_type, 
                              levels = c("no_carbohydrate", "rapid_digestible", "slow_digestible"))
  )

for (analyte_name in scfa_analytes) {
  case_analyte_data <- case_only_data %>% filter(analyte == analyte_name)
  
  # Model with time and carbohydrate effects, subject as random intercept
  time_carb_model <- lmer(concentration ~ timepoint_hr * carbohydrate_type + 
                          (1 | subject), 
                          data = case_analyte_data, 
                          REML = FALSE)
  
  # Model with time effect only, subject as random intercept
  time_only_model <- lmer(concentration ~ timepoint_hr + (1 | subject), 
                          data = case_analyte_data, 
                          REML = FALSE)
  
  # Model with carbohydrate effect only, subject as random intercept
  carb_only_model <- lmer(concentration ~ carbohydrate_type + (1 | subject), 
                          data = case_analyte_data, 
                          REML = FALSE)
  
  # Null model (intercept + random effect only)
  null_model <- lmer(concentration ~ 1 + (1 | subject), 
                     data = case_analyte_data, 
                     REML = FALSE)
  
  # Store results
  case_only_mixed_effects[[analyte_name]] <- list(
    time_carb_model = time_carb_model,
    time_only_model = time_only_model,
    carb_only_model = carb_only_model,
    null_model = null_model,
    anova_time_carb = anova(time_carb_model),
    summary_time_carb = summary(time_carb_model),
    # Model comparisons
    time_vs_null = anova(null_model, time_only_model),
    carb_vs_null = anova(null_model, carb_only_model),
    interaction_vs_main = anova(time_only_model, time_carb_model)
  )
}

# Save subject-level analyses to CSV files
write.csv(subject_summary, "results/subject_level_summary.csv", row.names = FALSE)
write.csv(within_subject_stats, "results/within_subject_changes.csv", row.names = FALSE)
write.csv(case_only_stats, "results/case_only_temporal_changes.csv", row.names = FALSE)
write.csv(case_only_pooled_stats, "results/case_only_pooled_temporal_changes.csv", row.names = FALSE)

cat("\n=== SUBJECT-LEVEL ANALYSES ===\n")
cat("\nSubject-level Summary Statistics:\n")
print(subject_summary)
cat("\nWithin-Subject Changes (0h to 48h):\n")
print(within_subject_stats)

cat("\n=== CASE-ONLY TEMPORAL ANALYSIS ===\n")
cat("\nCase Group Temporal Changes by Carbohydrate Type:\n")
print(case_only_stats)
cat("\nCase Group Temporal Changes (Pooled):\n")
print(case_only_pooled_stats)

cat("\n=== CASE-ONLY MIXED-EFFECTS MODEL RESULTS ===\n")
for (analyte_name in names(case_only_mixed_effects)) {
  cat("\n", analyte_name, " (Case Group Only):\n", sep = "")
  cat("ANOVA Results (Time × Carbohydrate interaction model):\n")
  print(case_only_mixed_effects[[analyte_name]]$anova_time_carb)
  
  cat("\nModel Comparisons:\n")
  cat("Time effect vs Null model:\n")
  print(case_only_mixed_effects[[analyte_name]]$time_vs_null)
  cat("Carbohydrate effect vs Null model:\n")  
  print(case_only_mixed_effects[[analyte_name]]$carb_vs_null)
  cat("Time × Carbohydrate interaction vs Time-only model:\n")
  print(case_only_mixed_effects[[analyte_name]]$interaction_vs_main)
  cat("\n" , rep("-", 50), "\n", sep = "")
}

# ===============================
# DELTA-FOCUSED ANALYSIS (48h - 0h changes)
# ===============================

# Calculate delta values for all subjects (more comprehensive than within_subject_changes)
delta_analysis <- sample_data_long %>%
  filter(timepoint_hr %in% c("0", "48")) %>%
  mutate(timepoint_hr = as.character(timepoint_hr)) %>%
  pivot_wider(names_from = timepoint_hr, 
              values_from = concentration, 
              names_prefix = "time_") %>%
  mutate(delta_48h_0h = time_48 - time_0) %>%
  filter(!is.na(delta_48h_0h))

# Statistical comparisons of delta values between groups
delta_group_stats <- delta_analysis %>%
  group_by(analyte) %>%
  t_test(delta_48h_0h ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Statistical comparisons of delta values between carbohydrate types
delta_carb_stats <- delta_analysis %>%
  group_by(analyte) %>%
  anova_test(delta_48h_0h ~ carbohydrate_type) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Post-hoc comparisons for carbohydrate effects on delta
delta_carb_posthoc <- delta_analysis %>%
  group_by(analyte) %>%
  pairwise_t_test(delta_48h_0h ~ carbohydrate_type, 
                  ref.group = "no_carbohydrate") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Two-way ANOVA for delta: Group × Carbohydrate
delta_interaction_stats <- delta_analysis %>%
  group_by(analyte) %>%
  anova_test(delta_48h_0h ~ group * carbohydrate_type) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Mixed-effects models for delta values
delta_mixed_effects <- list()

for (analyte_name in scfa_analytes) {
  delta_analyte_data <- delta_analysis %>% 
    filter(analyte == analyte_name) %>%
    mutate(
      group = factor(group, levels = c("control", "case")),
      carbohydrate_type = factor(carbohydrate_type, 
                                levels = c("no_carbohydrate", "rapid_digestible", "slow_digestible")),
      subject = factor(subject)
    )
  
  # Full interaction model
  full_model <- lmer(delta_48h_0h ~ group * carbohydrate_type + (1 | subject), 
                     data = delta_analyte_data, REML = FALSE)
  
  # Group only model
  group_model <- lmer(delta_48h_0h ~ group + (1 | subject), 
                      data = delta_analyte_data, REML = FALSE)
  
  # Carbohydrate only model
  carb_model <- lmer(delta_48h_0h ~ carbohydrate_type + (1 | subject), 
                     data = delta_analyte_data, REML = FALSE)
  
  # Null model
  null_model <- lmer(delta_48h_0h ~ 1 + (1 | subject), 
                     data = delta_analyte_data, REML = FALSE)
  
  delta_mixed_effects[[analyte_name]] <- list(
    full_model = full_model,
    group_model = group_model,
    carb_model = carb_model,
    null_model = null_model,
    anova_full = anova(full_model),
    summary_full = summary(full_model),
    group_vs_null = anova(null_model, group_model),
    carb_vs_null = anova(null_model, carb_model),
    interaction_vs_additive = anova(group_model, full_model)
  )
}

# Summary statistics for delta values by group and carbohydrate
delta_summary <- delta_analysis %>%
  group_by(group, carbohydrate_type, analyte) %>%
  summarise(
    n = n(),
    mean_delta = mean(delta_48h_0h, na.rm = TRUE),
    median_delta = median(delta_48h_0h, na.rm = TRUE),
    sd_delta = sd(delta_48h_0h, na.rm = TRUE),
    sem_delta = sd_delta / sqrt(n),
    q25_delta = quantile(delta_48h_0h, 0.25, na.rm = TRUE),
    q75_delta = quantile(delta_48h_0h, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Save delta analysis results
write.csv(delta_summary, "results/delta_summary_statistics.csv", row.names = FALSE)
write.csv(delta_group_stats, "results/delta_group_comparisons.csv", row.names = FALSE)
write.csv(delta_carb_stats, "results/delta_carbohydrate_effects.csv", row.names = FALSE)
write.csv(delta_carb_posthoc, "results/delta_carbohydrate_posthoc.csv", row.names = FALSE)

cat("\nCase-only analyses saved to results/ directory\n")
cat("\nSubject-level analyses saved to results/ directory\n")
cat("\n=== DELTA (48h-0h) ANALYSIS RESULTS ===\n")
cat("\nDelta Group Comparisons (Control vs Case Response Magnitude):\n")
print(delta_group_stats)
cat("\nDelta Carbohydrate Effects (Response Magnitude by Carb Type):\n")
print(delta_carb_stats)
cat("\nDelta Carbohydrate Post-hoc (vs No Carbohydrate):\n")
print(delta_carb_posthoc)

cat("\nDelta Group × Carbohydrate Interactions:\n")
delta_interaction_display <- as.data.frame(delta_interaction_stats) %>%
  select(analyte, Effect, p, p.adj, p.adj.signif)
print(delta_interaction_display)

cat("\nDelta Summary Statistics by Group and Carbohydrate:\n")
print(delta_summary)

cat("\n=== DELTA MIXED-EFFECTS MODEL RESULTS ===\n")
for (analyte_name in names(delta_mixed_effects)) {
  cat("\n", analyte_name, " (Delta Analysis):\n", sep = "")
  cat("ANOVA Results (Group × Carbohydrate model on delta values):\n")
  print(delta_mixed_effects[[analyte_name]]$anova_full)
  
  cat("\nModel Comparisons:\n")
  cat("Group effect vs Null model:\n")
  print(delta_mixed_effects[[analyte_name]]$group_vs_null)
  cat("Carbohydrate effect vs Null model:\n")  
  print(delta_mixed_effects[[analyte_name]]$carb_vs_null)
  cat("Group × Carbohydrate interaction:\n")
  print(delta_mixed_effects[[analyte_name]]$interaction_vs_additive)
  cat("\n" , rep("-", 50), "\n", sep = "")
}

cat("\nDelta (48h-0h) analyses saved to results/ directory\n")

cat("\n=== MIXED-EFFECTS MODEL RESULTS ===\n")
for (analyte_name in names(mixed_effects_results)) {
  cat("\n", analyte_name, ":\n", sep = "")
  print(mixed_effects_results[[analyte_name]]$anova_full)
}

# Print summary tables
cat("\nSummary by Group:\n")
print(summary_by_group)
cat("\nSummary by Carbohydrate Type:\n")
print(summary_by_carb)
cat("\nSummary by Timepoint:\n")
print(summary_by_time)
cat("\nCombined summary saved to results/combined_summary_statistics.csv\n")

# ===============================
# Visualization
# ===============================

# Define publication-quality theme
pub_theme <- theme_classic() +
  theme(
    text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# Custom color palette
group_colors <- c("control" = "#3498DB", "case" = "#E74C3C")
carb_colors <- c("no_carbohydrate" = "#95A5A6", 
                 "rapid_digestible" = "#F39C12", 
                 "slow_digestible" = "#27AE60")

# 1. Box plots by group faceted by carbohydrate type with statistical comparisons
plot_by_group <- sample_data_long %>%
  ggplot(aes(x = group, y = concentration, fill = group)) +
  geom_boxplot(alpha = 0.8) +
  facet_grid(carbohydrate_type ~ analyte, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("control", "case")),
                     label = "p.signif",
                     step_increase = 0.1) +
  labs(title = "SCFA Concentrations by Group and Carbohydrate Type",
       x = "Group",
       y = "Concentration (μM)",
       fill = "Group") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Box plots by carbohydrate type faceted by time with statistical comparisons
plot_by_carb <- sample_data_long %>%
  ggplot(aes(x = carbohydrate_type, y = concentration, fill = carbohydrate_type)) +
  geom_boxplot(alpha = 0.8) +
  facet_grid(timepoint_hr ~ analyte, scales = "free_y", 
             labeller = labeller(timepoint_hr = function(x) paste0(x, "h"))) +
  scale_fill_manual(values = carb_colors) +
  stat_compare_means(method = "t.test",
                     ref.group = "no_carbohydrate",
                     label = "p.signif",
                     step_increase = 0.1) +
  labs(title = "SCFA Concentrations by Carbohydrate Type and Time Point",
       x = "Carbohydrate Type",
       y = "Concentration (μM)",
       fill = "Carbohydrate Type") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# 3. Time series plots showing changes over time by carbohydrate type
plot_time_series <- sample_data_long %>%
  mutate(timepoint_hr = as.numeric(as.character(timepoint_hr))) %>%
  group_by(carbohydrate_type, timepoint_hr, analyte) %>%
  summarise(
    mean_conc = mean(concentration, na.rm = TRUE),
    sem = sd(concentration, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = timepoint_hr, y = mean_conc, color = carbohydrate_type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_conc - sem, ymax = mean_conc + sem), 
                width = 1, linewidth = 0.8) +
  facet_wrap(~ analyte, scales = "free_y") +
  scale_color_manual(values = carb_colors) +
  scale_x_continuous(breaks = c(0, 48)) +
  labs(title = "SCFA Concentrations Over Time by Carbohydrate Type",
       x = "Time (hours)",
       y = "Mean Concentration ± SEM (μM)",
       color = "Carbohydrate Type") +
  pub_theme +
  theme(legend.position = "bottom")

# 4. Heatmap showing mean concentrations across all conditions
heatmap_data <- sample_data_long %>%
  mutate(timepoint_hr = paste0(timepoint_hr, "h")) %>%
  group_by(group, carbohydrate_type, timepoint_hr, analyte) %>%
  summarise(mean_conc = mean(concentration, na.rm = TRUE), .groups = "drop") %>%
  unite("condition", group, carbohydrate_type, timepoint_hr, sep = "_")

plot_heatmap <- heatmap_data %>%
  ggplot(aes(x = condition, y = analyte, fill = mean_conc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(name = "Mean\nConcentration\n(μM)", option = "plasma") +
  labs(title = "SCFA Concentration Heatmap Across All Conditions",
       x = "Condition (Group_Carbohydrate_Time)",
       y = "SCFA Analyte") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "italic"))

# 5. Interaction plot showing group x carbohydrate effects over time
plot_interaction <- sample_data_long %>%
  mutate(timepoint_hr = as.numeric(as.character(timepoint_hr))) %>%
  group_by(group, carbohydrate_type, timepoint_hr, analyte) %>%
  summarise(mean_conc = mean(concentration, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = timepoint_hr, y = mean_conc, 
             color = carbohydrate_type, linetype = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, aes(shape = group)) +
  facet_wrap(~ analyte, scales = "free_y") +
  scale_color_manual(values = carb_colors) +
  scale_x_continuous(breaks = c(0, 48)) +
  scale_linetype_manual(values = c("control" = "solid", "case" = "dashed")) +
  scale_shape_manual(values = c("control" = 16, "case" = 17)) +
  labs(title = "SCFA Concentrations: Group × Carbohydrate × Time Interaction",
       x = "Time (hours)",
       y = "Mean Concentration (μM)",
       color = "Carbohydrate Type",
       linetype = "Group",
       shape = "Group") +
  pub_theme +
  theme(legend.position = "bottom")

# 6. Subject-level heatmap showing individual responses
subject_heatmap_data <- sample_data_long %>%
  mutate(
    condition_time = paste0(carbohydrate_type, "_", timepoint_hr, "h"),
    mean_conc = concentration  # Already averaged technical replicates
  ) %>%
  arrange(group, subject, analyte)

# Create ordered subject labels with control group first, then case group
subject_order <- subject_heatmap_data %>%
  select(subject, group) %>%
  distinct() %>%
  arrange(desc(group), subject) %>%  # control comes before case alphabetically
  mutate(
    subject_group = paste0(subject, " (", group, ")"),
    # Add spacing indicator for plotting
    group_spacing = ifelse(group == "control", "Control Group", "Case Group")
  )

# Create subject ordering with gap between groups
control_subjects <- subject_order %>% filter(group == "control") %>% arrange(subject)
case_subjects <- subject_order %>% filter(group == "case") %>% arrange(subject)

# Create spacer rows for the gap between groups
n_analytes <- length(unique(subject_heatmap_data$analyte))
n_conditions <- length(unique(subject_heatmap_data$condition_time))

# Create empty spacer data for gap
spacer_data <- expand_grid(
  analyte = unique(subject_heatmap_data$analyte),
  condition_time = unique(subject_heatmap_data$condition_time),
  subject_group = " ",
  mean_conc = NA_real_
)

# Create ordered factor levels with gap
all_subject_levels <- c(
  rev(control_subjects$subject_group),  # Control at top (reversed for ggplot)
  " ",           # Gap
  rev(case_subjects$subject_group)      # Case at bottom (reversed for ggplot)
)

# Add the ordered subject labels to the data and include spacer
subject_heatmap_data <- subject_heatmap_data %>%
  left_join(subject_order %>% select(subject, subject_group, group_spacing), by = "subject") %>%
  bind_rows(spacer_data) %>%
  mutate(
    subject_group = factor(subject_group, levels = all_subject_levels)
  )

plot_subject_heatmap <- subject_heatmap_data %>%
  ggplot(aes(x = condition_time, y = subject_group, fill = mean_conc)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ analyte, scales = "free", ncol = 2) +
  scale_fill_viridis_c(name = "Concentration\n(μM)", option = "plasma", 
                       trans = "sqrt", na.value = "white") +
  labs(title = "Individual Subject SCFA Responses Across All Conditions",
       subtitle = "Control group (top) and Case group (bottom) with visual separation",
       x = "Condition × Time Point",
       y = "Subject (Group)") +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(0.5, "lines")
  )

# 7. Individual subject trajectories by carbohydrate source
# Prepare data for trajectory plots
trajectory_data <- sample_data_long %>%
  mutate(
    timepoint_hr = as.numeric(as.character(timepoint_hr)),
    mean_conc = concentration  # Already averaged technical replicates
  )

# Calculate group means for overlay
group_means_trajectories <- trajectory_data %>%
  group_by(group, carbohydrate_type, timepoint_hr, analyte) %>%
  summarise(
    mean_conc = mean(mean_conc, na.rm = TRUE),
    sem = sd(mean_conc, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Create trajectory plots by carbohydrate source
plot_trajectories <- trajectory_data %>%
  ggplot(aes(x = timepoint_hr, y = mean_conc)) +
  # Individual subject lines (thin, semi-transparent)
  geom_line(aes(group = subject, color = group), 
            alpha = 0.4, linewidth = 0.5) +
  geom_point(aes(color = group), alpha = 0.5, size = 1) +
  # Group mean lines (thick, prominent)
  geom_line(data = group_means_trajectories, 
            aes(color = group), 
            linewidth = 2, alpha = 0.9) +
  geom_point(data = group_means_trajectories, 
             aes(color = group), 
             size = 3, alpha = 0.9) +
  # Error bars for group means
  geom_errorbar(data = group_means_trajectories,
                aes(ymin = mean_conc - sem, ymax = mean_conc + sem, color = group),
                width = 1, linewidth = 1, alpha = 0.7) +
  facet_grid(carbohydrate_type ~ analyte, scales = "free_y",
             labeller = labeller(carbohydrate_type = function(x) {
               str_replace_all(str_to_title(str_replace_all(x, "_", " ")), 
                              c("No" = "No", "Rapid" = "Rapid", "Slow" = "Slow"))
             })) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = c(0, 48)) +
  labs(title = "Individual Subject SCFA Trajectories by Carbohydrate Source",
       subtitle = "Thin lines: individual subjects, Thick lines: group means ± SEM",
       x = "Time (hours)",
       y = "Concentration (μM)",
       color = "Group") +
  pub_theme +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(0.3, "lines")
  )

# 8. Case-only temporal changes visualization
case_only_plot_data <- sample_data_long %>%
  filter(group == "case") %>%
  mutate(timepoint_hr = as.numeric(as.character(timepoint_hr))) %>%
  group_by(carbohydrate_type, timepoint_hr, analyte) %>%
  summarise(
    mean_conc = mean(concentration, na.rm = TRUE),
    sem = sd(concentration, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

plot_case_only <- case_only_plot_data %>%
  ggplot(aes(x = timepoint_hr, y = mean_conc, color = carbohydrate_type)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_conc - sem, ymax = mean_conc + sem), 
                width = 1, linewidth = 1) +
  facet_wrap(~ analyte, scales = "free_y") +
  scale_color_manual(values = carb_colors) +
  scale_x_continuous(breaks = c(0, 48)) +
  labs(title = "Case Group Only: SCFA Temporal Changes by Carbohydrate Type",
       subtitle = "Mean ± SEM concentrations from 0h to 48h",
       x = "Time (hours)",
       y = "Mean Concentration ± SEM (μM)",
       color = "Carbohydrate Type") +
  pub_theme +
  theme(legend.position = "bottom")

# Individual case subject trajectories with statistical overlay
case_individual_data <- sample_data_long %>%
  filter(group == "case") %>%
  mutate(
    timepoint_hr = as.numeric(as.character(timepoint_hr)),
    mean_conc = concentration  # Already averaged technical replicates
  )

plot_case_individual <- case_individual_data %>%
  ggplot(aes(x = timepoint_hr, y = mean_conc)) +
  # Individual subject lines
  geom_line(aes(group = subject), alpha = 0.6, linewidth = 0.8, color = "#E74C3C") +
  geom_point(alpha = 0.6, size = 2, color = "#E74C3C") +
  # Group mean overlay
  geom_line(data = case_only_plot_data, 
            aes(color = carbohydrate_type), 
            linewidth = 2, alpha = 0.9) +
  geom_point(data = case_only_plot_data, 
             aes(color = carbohydrate_type), 
             size = 4, alpha = 0.9) +
  facet_grid(carbohydrate_type ~ analyte, scales = "free_y",
             labeller = labeller(carbohydrate_type = function(x) {
               str_replace_all(str_to_title(str_replace_all(x, "_", " ")), 
                              c("No" = "No", "Rapid" = "Rapid", "Slow" = "Slow"))
             })) +
  scale_color_manual(values = carb_colors) +
  scale_x_continuous(breaks = c(0, 48)) +
  labs(title = "Case Group Individual Subject SCFA Trajectories",
       subtitle = "Thin lines: individual subjects, Thick lines: group means",
       x = "Time (hours)",
       y = "Concentration (μM)",
       color = "Carbohydrate Type") +
  pub_theme +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(0.3, "lines")
  )

# 8. Delta (48h-0h) focused visualizations
# Group comparison of response magnitudes
plot_delta_group <- delta_analysis %>%
  ggplot(aes(x = group, y = delta_48h_0h, fill = group)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~ analyte, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("control", "case")),
                     label = "p.signif",
                     step_increase = 0.1) +
  labs(title = "SCFA Response Magnitude (Δ48h-0h) by Group",
       subtitle = "Comparing metabolic response between control and case groups",
       x = "Group",
       y = "Delta Concentration (48h - 0h) μM",
       fill = "Group") +
  pub_theme +
  theme(legend.position = "bottom")

# Carbohydrate comparison of response magnitudes
plot_delta_carb <- delta_analysis %>%
  ggplot(aes(x = carbohydrate_type, y = delta_48h_0h, fill = carbohydrate_type)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~ analyte, scales = "free_y") +
  scale_fill_manual(values = carb_colors) +
  stat_compare_means(method = "t.test",
                     ref.group = "no_carbohydrate",
                     label = "p.signif",
                     step_increase = 0.1) +
  labs(title = "SCFA Response Magnitude (Δ48h-0h) by Carbohydrate Type",
       subtitle = "Comparing metabolic response magnitude across carbohydrate interventions",
       x = "Carbohydrate Type",
       y = "Delta Concentration (48h - 0h) μM",
       fill = "Carbohydrate Type") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Group × Carbohydrate interaction for delta
plot_delta_interaction <- delta_analysis %>%
  group_by(group, carbohydrate_type, analyte) %>%
  summarise(
    mean_delta = mean(delta_48h_0h, na.rm = TRUE),
    sem_delta = sd(delta_48h_0h, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = carbohydrate_type, y = mean_delta, fill = group)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_delta - sem_delta, ymax = mean_delta + sem_delta),
                position = position_dodge(width = 0.9), width = 0.2) +
  facet_wrap(~ analyte, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(title = "SCFA Response Magnitude: Group × Carbohydrate Interaction",
       subtitle = "Mean ± SEM change from baseline (48h - 0h)",
       x = "Carbohydrate Type",
       y = "Mean Delta Concentration ± SEM (μM)",
       fill = "Group") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Display plots
print(plot_by_group)
print(plot_by_carb)
print(plot_time_series)
print(plot_heatmap)
print(plot_interaction)
print(plot_subject_heatmap)
print(plot_trajectories)
print(plot_case_only)
print(plot_case_individual)
print(plot_delta_group)
print(plot_delta_carb)
print(plot_delta_interaction)

# Save plots
ggsave("plots/scfa_by_group.png", plot_by_group, width = 14, height = 10, dpi = 300)
ggsave("plots/scfa_by_carb.png", plot_by_carb, width = 14, height = 8, dpi = 300)
ggsave("plots/scfa_time_series.png", plot_time_series, width = 12, height = 8, dpi = 300)
ggsave("plots/scfa_heatmap.png", plot_heatmap, width = 14, height = 6, dpi = 300)
ggsave("plots/scfa_interaction.png", plot_interaction, width = 14, height = 10, dpi = 300)
ggsave("plots/scfa_subject_heatmap.png", plot_subject_heatmap, width = 16, height = 12, dpi = 300)
ggsave("plots/scfa_trajectories.png", plot_trajectories, width = 16, height = 12, dpi = 300)
ggsave("plots/scfa_case_only_temporal.png", plot_case_only, width = 12, height = 8, dpi = 300)
ggsave("plots/scfa_case_individual_trajectories.png", plot_case_individual, width = 16, height = 12, dpi = 300)
ggsave("plots/scfa_delta_group_comparison.png", plot_delta_group, width = 12, height = 8, dpi = 300)
ggsave("plots/scfa_delta_carbohydrate_effects.png", plot_delta_carb, width = 12, height = 8, dpi = 300)
ggsave("plots/scfa_delta_group_carb_interaction.png", plot_delta_interaction, width = 14, height = 10, dpi = 300)

