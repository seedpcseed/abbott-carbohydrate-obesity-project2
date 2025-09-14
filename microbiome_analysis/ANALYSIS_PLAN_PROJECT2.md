```markdown
# ABBOTT–LURIE COLLABORATION  
## Project 2: Microbiome Effects from *Ex Vivo* Carbohydrate Exposures  
### Complete Analysis Plan (for Statistician & Coding Assistant)

> **Purpose**  
> Quantify how adolescent gut microbiota (Healthy vs Obese) adapt *ex vivo* to rapid- vs slow-digestible carbohydrates over time, using 16S rDNA sequencing. Account for repeated measures and **biological replicates** per subject×condition×time.



## 1) Data & Variables

### 1.1 Inputs
- **Counts**: ASV/OTU table (rows = features, cols = samples).  
- **Taxonomy**: feature→(Kingdom…Species/ASV).  
- **Metadata**: one row per sample with at least:
  - `SampleID` (rownames)  
  - `subjectID` (e.g., `84B`)  
  - `group` (levels: `Healthy`, `Obese`)  
  - `carbohydrate` (levels: `No Added Carb`, `RDC`, `SDC`)  
  - `timepoint` (levels: `0H`, `48H`, `72H`)  
  - `replicate` (integer 1/2; biological replicate indicator)  
  - optional: `raw_seq_file`, run/batch if needed for QC.

**Sanity**: `colnames(counts) == rownames(metadata)` after matching.

### 1.2 Factor Coding (set explicit reference levels)
- `group`: `Healthy` (ref) vs `Obese`  
- `carbohydrate`: `No Added Carb` (ref), `RDC`, `SDC`  
- `timepoint`: `0H` (ref), `48H`, `72H`

**R snippet**:
```r
md <- metadata %>% mutate(
  group        = factor(group,        levels = c("Healthy","Obese")),
  carbohydrate = factor(carbohydrate, levels = c("No Added Carb","RDC","SDC")),
  timepoint    = factor(timepoint,    levels = c("0H","48H","72H")),
  subjectID    = factor(subjectID),
  replicate    = factor(replicate)
)
````

---

## 2) Biological Replicates & Repeated Measures

### 2.1 Structure

* **Repeated measures**: same `subjectID` across `timepoint` and `carbohydrate`.
* **Biological replicates**: multiple wells per **subject×carbohydrate×timepoint** (e.g., R1/R2, S1/S2, NC).

### 2.2 Modeling Approach (**Primary**)

* Treat **subject** as the primary clustering unit.
* Treat **biological replicates** as random intercepts **nested** within subject×condition×time.

**Keys**:

```r
md <- md %>%
  mutate(cell = interaction(subjectID, carbohydrate, timepoint, drop = TRUE),
         replicate_id = interaction(subjectID, carbohydrate, timepoint, replicate, drop = TRUE))
```

* Random-effects formulas (examples):

  * Alpha: `(1|subjectID) + (1|subjectID:cell)`, or `(1|subjectID) + (1|subjectID:replicate_id)`
  * DA (ANCOM-BC2): `rand_formula = (1|subjectID) + (1|subjectID:cell)`

### 2.3 Sensitivity / Visualization

* **Collapse replicates** to *subject×carbohydrate×timepoint* by **summing counts** (then renormalize) for:

  * Cleaner ordinations/heatmaps.
  * Sensitivity PERMANOVA and DA (documented).

---

## 3) Quality Control

### 3.1 Sequencing Depth & Filtering

* Remove samples with read depth `< Q10` of library sizes or `<1,000` reads (pre-specify).
* Filter features with total counts `< 10` and prevalence `< 5%` (pre-analysis visualization can be looser).

### 3.2 Replicate Concordance (before any collapse)

* Within each `cell` (subject×carb×time):

  * **Bray–Curtis** pairwise distances among replicates; flag if `> 95th` percentile of within-cell distances.
  * **Spearman ρ** of genus-level relative abundance; flag if `< 0.6`.
* **Decision**:

  * If one replicate fails QC → **exclude that replicate**.
  * If both pass → **keep both** (primary mixed-effects analysis).
  * Document all exclusions in `outputs/tables/qc_exclusions.csv`.


## 4) Phyloseq Assembly & Normalization

* Build `phyloseq` object (counts, taxonomy, metadata).
* For **beta**/ordination: use **relative abundance** (`transform_sample_counts(x/sum(x))`).
* For **DA**: use methods’ native handling (ANCOM-BC2/CLR); avoid rarefaction for inference.

---

## 5) Alpha Diversity

### 5.1 Metrics

* Observed richness, Shannon, Simpson, InvSimpson, Evenness (Shannon/log(Observed)).

### 5.2 Model (per metric)

* **Fixed**: `group * carbohydrate * timepoint`
* **Random**: `(1|subjectID) + (1|subjectID:cell)`
* **Diagnostics**: residual plots; if heteroskedastic, use `nlme::lme` with `varIdent(~1|timepoint)`.

**Example**:

```r
fit <- lme4::lmer(Shannon ~ group*carbohydrate*timepoint +
                    (1|subjectID) + (1|subjectID:cell),
                  data = alpha_df)
anova_tbl <- car::Anova(fit, type="III")
emm <- emmeans::emmeans(fit, ~ carbohydrate | timepoint * group)
pairs_emm <- summary(pairs(emm), adjust="fdr")
```

**Report**: effects (β), 95% CI, p, FDR-adjusted for post-hoc families.

---

## 6) Beta Diversity

### 6.1 Distances & Ordination

* **Primary**: Bray–Curtis and Jaccard on relative abundance.
* **Sensitivity**: Weighted/Unweighted UniFrac if a tree is available.
* Ordination: PCoA; plot by `carbohydrate×timepoint`, facet by `group`.

### 6.2 PERMANOVA (primary inference)

* **Formula**: `distance ~ group * carbohydrate * timepoint`
* **Strata**: permutations **constrained within `subjectID`** (respect repeated measures).
* **Replicates**: keep both (post-QC); permutations remain within subject.

```r
set.seed(42)
permanova <- vegan::adonis2(dist ~ group*carbohydrate*timepoint,
                            data   = md,
                            strata = md$subjectID,
                            permutations = 999)
```

* **Pairwise** only if omnibus significant; FDR-adjust p’s within families.
* **Effect size**: report R² and CI (via permutation or bootstrap).

### 6.3 Sensitivity (collapsed replicates)

* Aggregate by `cell` (sum counts), re-run PERMANOVA; convergence of results supports robustness.

---

## 7) Differential Abundance (DA)

### 7.1 Primary: **ANCOM-BC2** with Random Effects

```r
res_ab <- ANCOMBC::ancombc2(
  data         = ps,  # phyloseq with all replicates kept
  formula      = "group * carbohydrate * timepoint",
  rand_formula = "(1|subjectID) + (1|subjectID:cell)",
  global       = TRUE,
  p_adj_method = "fdr",
  prv_cut = 0.10, lib_cut = 1000, alpha = 0.05, verbose = FALSE
)
```

* **Global tests** per multi-level factor; if significant, follow with **simple contrasts**:

  * `RDC vs No Added Carb`, `SDC vs No Added Carb`, `RDC vs SDC` at each timepoint and by group where relevant.
  * `48H vs 0H`, `72H vs 0H` within each carb & group.
* **Outputs**: effect (log-fold change), SE, 95% CI, q-value, prevalence.

### 7.2 Secondary (Sensitivity)

* **MaAsLin2** with `random_effects=c("subjectID","cell")`, normalization=`CLR` or `TSS`.
* **DESeq2** blocking by subject (fixed effect): `~ subjectID + group*carbohydrate*timepoint` (note: no random effects).
* **Corncob** only after **collapsing** replicates (or as exploratory), since it lacks random effects.

### 7.3 Multiple Testing

* BH-FDR within each coherent family (e.g., all taxa for one contrast).
* Report both raw and adjusted p; emphasize **q<0.05**.

---

## 8) Targeted Contrasts (Pre-Specified)

1. **Time Adaptation** (within each carb & group):

   * `48H vs 0H`, `72H vs 0H`
2. **Carbohydrate Type** (at each timepoint & group):

   * `RDC vs SDC`, `RDC vs No Added Carb`, `SDC vs No Added Carb`
3. **Group Differences**:

   * Obese vs Healthy at baseline (0H) and trajectory differences (group×time×carb interactions; report simple slopes/contrasts).
4. **Global Interactions**:

   * If `group*carbohydrate*timepoint` significant, present stratified effects.

---

## 9) Visualization Standards

* **Alpha**: violin/box + jitter; color by `carbohydrate`, facet by `group`; annotate EMM differences.
* **Beta**: PCoA scatter; hulls/ellipses per `carbohydrate×timepoint`; facet by `group`; caption includes PERMANOVA R²/p.
* **DA**:

  * Volcano/Manhattan: x = effect size, y = –log10(q); highlight `q<0.05`.
  * Heatmap of significant taxa (prevalence ≥20%, mean rel. abund ≥0.1%); columns ordered by `group→carb→time`, rows clustered.

---

## 10) Inline, Dynamic Interpretation (Rmd)

* No hard-coded numbers. Use helper functions to render sentences from current model outputs.

**Examples**:

```r
interp_alpha <- function(fit, metric){
  an <- car::Anova(fit, type="III")
  p <- an["timepoint","Pr(>Chisq)"]
  if(!is.na(p) && p < 0.05) {
    glue::glue("{metric} changed over time (Type III p = {signif(p,3)}).")
  } else glue::glue("No significant time effect on {metric}.")
}
```

```r
interp_permanova <- function(perma){
  r2 <- perma$R2[rownames(perma)=="group:carbohydrate:timepoint"]
  p  <- perma$`Pr(>F)`[rownames(perma)=="group:carbohydrate:timepoint"]
  if(length(p)==1 && !is.na(p) && p<0.05)
    glue::glue("Community composition differed by the group×carbohydrate×timepoint interaction (R² = {round(r2,3)}, p = {signif(p,3)}).")
  else "No significant three-way interaction; see main effects."
}
```

---

## 11) Output Artifacts

* **Tables (CSV + nicely formatted HTML):**

  * `alpha_mixed_models_[metric].csv` (fixed effects, random effects variance).
  * `alpha_emmeans_posthoc.csv` (FDR-adjusted).
  * `permanova_summary.csv` (effects, R², p).
  * `da_ancombc2_[contrast].csv` (taxon, effect, CI, q).
  * `qc_summary.csv`, `qc_exclusions.csv`, `cell_counts.csv`.

* **Figures (PNG + PDF)**:

  * `alpha_[metric]_[facet].(png|pdf)`
  * `pcoa_[distance]_[facet].(png|pdf)`
  * `da_volcano_[contrast].(png|pdf)`
  * `da_heatmap_sig_taxa.(png|pdf)`
  * `rarefaction_curves.(png|pdf)` (QC)

* **R objects (RDS)**:

  * `ps_raw.rds`, `ps_filtered.rds`, `ps_rel.rds`, `ps_collapsed.rds`
  * `fit_alpha_[metric].rds`, `permanova_[distance].rds`, `ancombc2_results.rds`

---

## 12) Power & Sensitivity (If Feasible)

* Post-hoc power for alpha (variance from mixed models) and for PERMANOVA (via simulation/permutation).
* Note expected effect sizes from prior rat work; emphasize confidence intervals in primary reporting.

---

## 13) R Markdown Skeleton (Drop-In)

```yaml
---
title: "Ex Vivo Carbohydrate–Microbiome Adaptation: 16S Analysis"
output: html_document
params:
  min_reads: 1000
  feature_min_total: 10
  prevalence_min: 0.05
---
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=TRUE, message=FALSE, warning=FALSE, fig.width=9, fig.height=5, dpi=300)
set.seed(20250101)
library(tidyverse); library(phyloseq); library(vegan); library(microbiome)
library(lme4); library(nlme); library(car); library(emmeans); library(broom); library(broom.mixed)
library(ANCOMBC); library(patchwork); library(glue); library(kableExtra); library(viridis)
```

```{r data-import}
counts <- read.table("data_raw/otu_table.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
tax    <- read.table("data_raw/taxonomy_table.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
md     <- read.table("data_raw/metadata.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE) %>% as_tibble(rownames="SampleID")

# Filter
libsize <- colSums(counts)
keep_samp <- names(libsize)[libsize >= params$min_reads]
counts <- counts[, intersect(keep_samp, md$SampleID)]
md     <- md %>% filter(SampleID %in% colnames(counts))

# Factor levels
md <- md %>% mutate(
  group        = factor(group,        levels=c("Healthy","Obese")),
  carbohydrate = factor(carbohydrate, levels=c("No Added Carb","RDC","SDC")),
  timepoint    = factor(timepoint,    levels=c("0H","48H","72H")),
  subjectID    = factor(subjectID),
  replicate    = factor(replicate)
) %>% mutate(
  cell         = interaction(subjectID, carbohydrate, timepoint, drop=TRUE),
  replicate_id = interaction(subjectID, carbohydrate, timepoint, replicate, drop=TRUE)
)

# Feature filter
keep_feat <- rowSums(counts) >= params$feature_min_total
counts <- counts[keep_feat, ]

# Assemble phyloseq
OTU  <- otu_table(as.matrix(counts), taxa_are_rows=TRUE)
TAX  <- tax_table(as.matrix(tax[rownames(counts), , drop=FALSE]))
META <- sample_data(as.data.frame(md) %>% column_to_rownames("SampleID"))
ps   <- phyloseq(OTU, TAX, META)
```

```{r qc-replicates}
# Genus-level for QC
ps_genus <- tax_glom(ps, "Genus")
ps_rel   <- transform_sample_counts(ps_genus, function(x) x/sum(x))
md_rel   <- as.data.frame(sample_data(ps_rel)); md_rel$SampleID <- rownames(md_rel)

# Within-cell Bray distances among replicates
dist_bray <- phyloseq::distance(ps_rel, "bray")
get_within <- function(cell_id){
  sids <- md_rel$SampleID[md_rel$cell==cell_id]
  if(length(sids) < 2) return(NA_real_)
  as.matrix(dist_bray)[sids, sids] %>% .[upper.tri(., diag=FALSE)] %>% mean()
}
wd <- sapply(unique(md_rel$cell), get_within)
thresh <- quantile(wd, 0.95, na.rm=TRUE)

# Flag samples in outlier cells (for review)
qc_cells_flag <- names(wd)[wd > thresh]
write.csv(data.frame(cell=names(wd), mean_within_bray=wd),
          "outputs/tables/qc_within_cell_bray.csv", row.names=FALSE)
```

```{r alpha}
alpha_metrics <- c("Observed","Shannon","Simpson","InvSimpson")
alpha_df <- estimate_richness(ps, measures=alpha_metrics) %>%
  rownames_to_column("SampleID") %>%
  left_join(as.data.frame(sample_data(ps)) %>% rownames_to_column("SampleID"),
            by="SampleID")

fits <- map(alpha_metrics, \(m){
  fm <- as.formula(paste0(m," ~ group*carbohydrate*timepoint + (1|subjectID) + (1|subjectID:cell)"))
  lme4::lmer(fm, data=alpha_df)
}) %>% set_names(alpha_metrics)

anova_tbl <- imap_dfr(fits, \(fit, m){
  aa <- car::Anova(fit, type="III") %>% broom::tidy() %>% mutate(metric=m)
})
write.csv(anova_tbl, "outputs/tables/alpha_anova_typeIII.csv", row.names=FALSE)
```

```{r beta}
ps_rel_all <- transform_sample_counts(ps, function(x) x/sum(x))
bray <- phyloseq::distance(ps_rel_all, "bray")
md_b  <- as.data.frame(sample_data(ps_rel_all))
set.seed(42)
permanova <- vegan::adonis2(bray ~ group*carbohydrate*timepoint,
                            data=md_b, strata=md_b$subjectID, permutations=999)
capture.output(permanova, file="outputs/tables/permanova_bray.txt")

ord <- ordinate(ps_rel_all, method="PCoA", distance="bray")
p_ord <- plot_ordination(ps_rel_all, ord, color="carbohydrate", shape="timepoint") +
  geom_point(size=2, alpha=0.8) + facet_wrap(~group) + theme_bw()
ggsave("outputs/figures/pcoa_bray.png", p_ord, width=8, height=5, dpi=300)
```

```{r da_ancombc2}
res_ab <- ANCOMBC::ancombc2(
  data         = ps,
  formula      = "group * carbohydrate * timepoint",
  rand_formula = "(1|subjectID) + (1|subjectID:cell)",
  global       = TRUE,
  p_adj_method = "fdr",
  prv_cut = 0.10, lib_cut = 1000, alpha = 0.05, verbose = FALSE
)
saveRDS(res_ab, "outputs/da/ancombc2_results.rds")

# Example: extract significant taxa for carbohydrate main effect (global or simple)
res_df <- res_ab$res %>% as_tibble()
write.csv(res_df, "outputs/tables/ancombc2_all_results.csv", row.names=FALSE)
```

```{r dynamic-interpretation}
interp_permanova <- function(perma){
  if(!"group:carbohydrate:timepoint" %in% rownames(perma)) return("PERMANOVA summary unavailable.")
  r2 <- perma$R2["group:carbohydrate:timepoint"]; p <- perma$`Pr(>F)`["group:carbohydrate:timepoint"]
  if(!is.na(p) && p < 0.05) glue::glue("Community composition differed by the group×carbohydrate×timepoint interaction (R²={round(r2,3)}, p={signif(p,3)}).")
  else "No significant three-way interaction; see main effects."
}

cat(interp_permanova(permanova))
```

```{r session}
sessionInfo()
```

---

## 14) Statistical Reporting Conventions

* **Alpha**: Type III tests for fixed effects; EMMs with FDR-adjusted contrasts; report β, CI, p, q.
* **Beta**: PERMANOVA (strata = subjectID), report R², p; pairwise only if omnibus significant; FDR across pairs.
* **DA**: ANCOM-BC2 effects with 95% CI and q; report direction and magnitude; include prevalence and mean rel. abundance.
* **Plots**: include N per panel; annotate significant contrasts (symbols or text with q).

---

## 15) Risk & Mitigation

* **Unbalanced cells / missing timepoints**: mixed models handle; report cell counts.
* **Low-depth samples**: filtered by min reads; sensitivity with/without borderlines.
* **Batch effects**: explore; include as random/fixed covariate if needed (sensitivity).
* **Replicate discordance**: pre-specified QC rules and exclusions table.

---

## 16) Deliverables

1. **HTML report** with dynamic text, figures, and linked CSV tables.
2. **All tables** (CSV) and **figures** (PNG/PDF) under `outputs/`.
3. **R objects** (RDS) for ps and model fits.
4. **QC appendix** documenting replicate checks and exclusions.
5. **README** detailing run instructions and package versions (`renv.lock`).

---

## 17) What to Hand Back If Data Change

* Re-knit the Rmd; all numbers auto-update (no hard-coding).
* Ensure factor levels remain consistent with specified reference levels.
* Re-run QC thresholds; update exclusions and note deltas in the report.

---

```
```
