# this script is to test the different parameters of the dada2-min-abundance

library(tidyverse)
library(readr)

# read in the asv tables
dada1 <- read_tsv("dada2_asv_1/asv_table.tsv") %>% 
select(-Sequence, -Length, -Total_Abundance, -Relative_Abundance, -ASV_ID)
dada2 <- read_tsv("dada2_asv_2/asv_table.tsv") %>% select(-Sequence, -Length, -Total_Abundance, -Relative_Abundance, -ASV_ID)
dada3 <- read_tsv("dada2_asv_3/asv_table.tsv") %>% select(-Sequence, -Length, -Total_Abundance, -Relative_Abundance, -ASV_ID)
dada4 <- read_tsv("dada2_asv_4/asv_table.tsv") %>% select(-Sequence, -Length, -Total_Abundance, -Relative_Abundance, -ASV_ID)
dada5 <- read_tsv("dada2_asv_5/asv_table.tsv") %>% select(-Sequence, -Length, -Total_Abundance, -Relative_Abundance, -ASV_ID)
dada6 <- read_tsv("dada2_asv_6/asv_table.tsv") %>% select(-Sequence, -Length, -Total_Abundance, -Relative_Abundance, -ASV_ID)
dada7 <- read_tsv("dada2_asv_7/asv_table.tsv") %>% select(-Sequence, -Length, -Total_Abundance, -Relative_Abundance, -ASV_ID)

# read in the taxonomy tables
tax1 <- read_tsv("dada2_asv_1/strobelca_taxonomy.tsv") %>% select(-Sequence, -Confidence, -Method, -Identity, -Coverage, -Reference_ID, -ASV_ID)
tax2 <- read_tsv("dada2_asv_2/strobelca_taxonomy.tsv") %>% select(-Sequence, -Confidence, -Method, -Identity, -Coverage, -Reference_ID, -ASV_ID)
tax3 <- read_tsv("dada2_asv_3/strobelca_taxonomy.tsv") %>% select(-Sequence, -Confidence, -Method, -Identity, -Coverage, -Reference_ID, -ASV_ID)
tax4 <- read_tsv("dada2_asv_4/strobelca_taxonomy.tsv") %>% select(-Sequence, -Confidence, -Method, -Identity, -Coverage, -Reference_ID, -ASV_ID)
tax5 <- read_tsv("dada2_asv_5/strobelca_taxonomy.tsv") %>% select(-Sequence, -Confidence, -Method, -Identity, -Coverage, -Reference_ID, -ASV_ID)
tax6 <- read_tsv("dada2_asv_6/strobelca_taxonomy.tsv") %>% select(-Sequence, -Confidence, -Method, -Identity, -Coverage, -Reference_ID, -ASV_ID)
tax7 <- read_tsv("dada2_asv_7/strobelca_taxonomy.tsv") %>% select(-Sequence, -Confidence, -Method, -Identity, -Coverage, -Reference_ID, -ASV_ID)

df <- data.frame(dada1 = nrow(dada1), dada2 = nrow(dada2), dada3 = nrow(dada3), dada4 = nrow(dada4), dada5 = nrow(dada5), dada6 = nrow(dada6), dada7 = nrow(dada7)) %>% t() %>% as.data.frame()
colnames(df) <- c("counts")
df <- df %>% rownames_to_column("dada2-min-abundance")

# now plot the number of ASVs for each dada2-min-abundance
ggplot(data = df, aes(x = `dada2-min-abundance`, y = counts)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "DADA2 min abundance", y = "Number of ASVs")

df_colsums <- data.frame(dada1 = colSums(dada1), dada2 = colSums(dada2), dada3 = colSums(dada3), dada4 = colSums(dada4), dada5 = colSums(dada5), dada6 = colSums(dada6), dada7 = colSums(dada7)) %>% t() %>% as.data.frame()
df_colsums <- df_colsums %>% rownames_to_column("dada2-min-abundance")


# Reshape data for better plotting
df_colsums_long <- df_colsums %>%
  pivot_longer(cols = -`dada2-min-abundance`, 
               names_to = "sample", 
               values_to = "total_counts")

# 1. Boxplot showing distribution of total counts across samples for each dada-min-abundance setting
p1 <- ggplot(data = df_colsums_long, aes(x = `dada2-min-abundance`, y = total_counts)) +
  geom_boxplot(aes(fill = `dada2-min-abundance`), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  theme_minimal() +
  labs(x = "DADA2 min abundance parameter", 
       y = "Total counts per sample",
       title = "Distribution of Total Sequencing Depth by DADA2 min-abundance Parameter") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# 2. Violin plot for better distribution visualization
p2 <- ggplot(data = df_colsums_long, aes(x = `dada2-min-abundance`, y = total_counts)) +
  geom_violin(aes(fill = `dada2-min-abundance`), alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.8) +
  theme_minimal() +
  labs(x = "DADA2 min abundance parameter", 
       y = "Total counts per sample",
       title = "Distribution of Total Sequencing Depth (Violin Plot)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# 3. Create a combined plot showing both ASV counts and total sequencing depth
df_combined <- df %>%
  left_join(df_colsums_long %>% 
              group_by(`dada2-min-abundance`) %>% 
              summarise(mean_counts = mean(total_counts),
                        median_counts = median(total_counts),
                        .groups = 'drop'),
            by = "dada2-min-abundance")

# Reshape for dual y-axis plotting
df_plot <- df_combined %>%
  select(`dada2-min-abundance`, counts, mean_counts) %>%
  pivot_longer(cols = c(counts, mean_counts), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric = case_when(
    metric == "counts" ~ "Number of ASVs",
    metric == "mean_counts" ~ "Mean Total Counts"
  ))

p3 <- ggplot(data = df_plot, aes(x = `dada2-min-abundance`, y = value, fill = metric)) +
  geom_col(position = "dodge", alpha = 0.8) +
  theme_minimal() +
  labs(x = "DADA2 min abundance parameter", 
       y = "Value",
       title = "Comparison of ASV Counts vs Mean Total Sequencing Depth",
       fill = "Metric") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Number of ASVs" = "steelblue", 
                               "Mean Total Counts" = "darkorange"))

print(p3)

# 4. Create summary statistics table
summary_stats <- df_colsums_long %>%
  group_by(`dada2-min-abundance`) %>%
  summarise(
    n_samples = n(),
    mean_counts = round(mean(total_counts), 0),
    median_counts = round(median(total_counts), 0),
    sd_counts = round(sd(total_counts), 0),
    min_counts = min(total_counts),
    max_counts = max(total_counts),
    .groups = 'drop'
  ) %>%
  left_join(df %>% select(`dada2-min-abundance`, counts), by = "dada2-min-abundance") %>%
  rename(n_asvs = counts)

print("Summary Statistics:")
print(summary_stats)

# shannon and observed species alpha diversity comparisons between the different dada2-min-abundance parameters
library(phyloseq)

dada1_ps <- phyloseq(otu_table(dada1, taxa_are_rows = TRUE), tax_table(tax1))
dada2_ps <- phyloseq(otu_table(dada2, taxa_are_rows = TRUE), tax_table(tax2))
dada3_ps <- phyloseq(otu_table(dada3, taxa_are_rows = TRUE), tax_table(tax3))
dada4_ps <- phyloseq(otu_table(dada4, taxa_are_rows = TRUE), tax_table(tax4))
dada5_ps <- phyloseq(otu_table(dada5, taxa_are_rows = TRUE), tax_table(tax5))
dada6_ps <- phyloseq(otu_table(dada6, taxa_are_rows = TRUE), tax_table(tax6))
dada7_ps <- phyloseq(otu_table(dada7, taxa_are_rows = TRUE), tax_table(tax7))

dada1_alpha <- estimate_richness(dada1_ps, measures = c("Shannon", "Observed")) %>% 
    mutate(min_abundance = 1)
dada2_alpha <- estimate_richness(dada2_ps, measures = c("Shannon", "Observed")) %>% 
    mutate(min_abundance = 2)
dada3_alpha <- estimate_richness(dada3_ps, measures = c("Shannon", "Observed")) %>% 
    mutate(min_abundance = 3)   
dada4_alpha <- estimate_richness(dada4_ps, measures = c("Shannon", "Observed")) %>% 
    mutate(min_abundance = 4)
dada5_alpha <- estimate_richness(dada5_ps, measures = c("Shannon", "Observed")) %>% 
    mutate(min_abundance = 5)
dada6_alpha <- estimate_richness(dada6_ps, measures = c("Shannon", "Observed")) %>% 
    mutate(min_abundance = 6)
dada7_alpha <- estimate_richness(dada7_ps, measures = c("Shannon", "Observed")) %>% 
    mutate(min_abundance = 7)

df_alpha_shannon <- bind_rows(dada1_alpha %>% select(-Observed), dada2_alpha %>% select(-Observed), dada3_alpha %>% select(-Observed), dada4_alpha %>% select(-Observed), dada5_alpha %>% select(-Observed), dada6_alpha %>% select(-Observed), dada7_alpha %>% select(-Observed)) %>%
    mutate(min_abundance = as.factor(min_abundance))
    
df_alpha_observed <- bind_rows(dada1_alpha %>% select(-Shannon), dada2_alpha %>% select(-Shannon), dada3_alpha %>% select(-Shannon), dada4_alpha %>% select(-Shannon), dada5_alpha %>% select(-Shannon), dada6_alpha %>% select(-Shannon), dada7_alpha %>% select(-Shannon)) %>%
    mutate(min_abundance = as.factor(min_abundance))

shannon_plot <-ggplot(df_alpha_shannon, aes(min_abundance, Shannon, fill = min_abundance)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "DADA2 min abundance parameter", y = "Shannon Diversity")

observed_plot <- ggplot(df_alpha_observed, aes(min_abundance, Observed, fill = min_abundance)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "DADA2 min abundance parameter", y = "Observed Species")

cowplot::plot_grid(shannon_plot, observed_plot, nrow = 1)