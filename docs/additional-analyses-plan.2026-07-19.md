# Additional analyses plan: qPCR-scaled microbiome, SCFA integration, and modern methods

**Project:** Abbott–Lurie ZC07, Project 2  
**Date:** July 19, 2026  
**Status:** Partially executed 2026-07-19 — Priority A/B core analyses completed under `integrated/results/additional_analyses_2026-07-19/`; MaAsLin 3 / ALDEx3 and facility qPCR QC remain blocked.  
**Prerequisite:** Retain the locked donor → A/B aliquot → culture well → time hierarchy from [`unit-hierarchy-reanalysis-plan.2026-07-19.md`](unit-hierarchy-reanalysis-plan.2026-07-19.md).  
**Execution log:** [`integrated/results/additional_analyses_2026-07-19/README.md`](../integrated/results/additional_analyses_2026-07-19/README.md)

## 1. Executive decision

Add a focused qPCR-scaled *Fusicatenibacter*–SCFA analysis, but do not run another simple sample-level correlation or pooled responder test.

The highest-value new analysis is:

> Is the SDC-associated 0-to-48-hour change in qPCR-scaled *Fusicatenibacter* abundance associated with the SDC-associated change in butyrate across independent donors?

The claim-facing analysis should use one value per donor (N=16), with exact well-level matching before replicate aggregation. An information-preserving mixed model using aliquot-level records can be a sensitivity analysis. The primary result must be a continuous association with an effect estimate and confidence interval, not a median responder comparison.

This is an exploratory association, not evidence that *Fusicatenibacter* directly consumes SDC or produces butyrate.

## 2. What the review found

### 2.1 The qPCR data support a quantitative analysis

The available files contain:

- 297 qPCR rows in `absolute.abundance.csv`;
- 292 qPCR-matched study profiles in the current phyloseq workflow;
- qPCR Ct values ranging from 9.46 to 38.62;
- reported `genome_copies_per_ul*` values ranging from 2 to approximately 6.60×10^9;
- *Fusicatenibacter* detected in 257 of 293 parseable study profiles (87.7%);
- *Fusicatenibacter* detected in 56 of 64 48-hour SDC profiles, representing 15 of 16 donors.

These counts support an abundance analysis. They do not support treating the 56 positive SDC profiles as 56 independent observations; the independent N remains 16 donors.

### 2.2 “Absolute Fusicatenibacter” is currently an indirect estimate

The current workflow calculates:

```text
qPCR-scaled Fusicatenibacter
  = Fusicatenibacter 16S relative abundance
  × total bacterial qPCR load
```

The dataset does not appear to contain a *Fusicatenibacter*-specific qPCR assay. Therefore:

- use **“qPCR-scaled *Fusicatenibacter* abundance”** or **“estimated *Fusicatenibacter* 16S target equivalents per µL”**;
- do not call the value a direct genus-specific qPCR measurement;
- do not call it cell count without a defensible 16S gene-copy correction;
- reserve direct genus/species quantification for a future targeted qPCR or ddPCR validation assay.

The source file reports both `gene_copies_per_ul` and `genome_copies_per_ul*`. The latter is exactly one quarter of the former in inspected rows, implying an assumed four 16S copies per genome. That assumption and the asterisk footnote must be obtained from the facility before manuscript use. Until then, use measured `gene_copies_per_ul` for scaling and describe the result as target equivalents rather than genomes or cells.

### 2.3 The current correlation is not claim-facing

The existing `fusicatenibacter-scfa-correlations` chunk:

1. uses relative rather than qPCR-scaled abundance;
2. restricts to 48-hour SDC samples;
3. runs ordinary `lm()` and Spearman tests on sample rows;
4. does not account for donor, aliquot, or well nesting;
5. reports row count as N even though the independent N is 16.

More importantly, the integration step first averages SCFA over wells within aliquot × carbohydrate × time, then joins that same average back onto individual microbiome well rows. That duplicates the SCFA outcome and prevents a true matched-well association.

The existing pooled absolute-abundance responder Wilcoxon analysis has a related problem: multiple records per donor are treated as independent and a continuous SCFA response is dichotomized.

Both analyses should be retired from claim-facing use when the replacement analysis is implemented.

### 2.4 The biological interpretation must remain open

The type-strain description of *Fusicatenibacter saccharivorans* reported lactate, formate, acetate, and succinate as glucose fermentation products, not butyrate. A positive association with butyrate could therefore reflect cross-feeding, co-expansion with another organism, total biomass, or shared response to substrate rather than direct butyrate production.

The analysis must compare butyrate with acetate, propionate, total SCFA, total bacterial load, and relative abundance before using pathway-specific language.

## 3. Analysis hierarchy

### Priority A: required before adding a manuscript claim

1. Complete adjusted 0-to-48-hour carbohydrate-change contrasts for acetate, propionate, and butyrate, including every pairwise carbohydrate contrast.
2. Replace claim-facing alpha-diversity Wilcoxon tests with donor-aware mixed models or paired donor-cell analyses.
3. Complete term-level PERMANOVA, donor-restricted permutations, and multivariate-dispersion testing.
4. Reconcile the 20 planned versus 16 analyzed donors and all missing 48-hour wells.
5. qPCR assay, unit, LOD/LOQ, and sample-match audit.
6. Exact SCFA–16S–qPCR join on donor, aliquot, carbohydrate, well, and time.
7. Total bacterial-load trajectory model.
8. Continuous qPCR-scaled *Fusicatenibacter*–butyrate analysis.
9. Donor-bootstrap and leave-one-donor-out influence analysis.
10. Relative-versus-qPCR-scaled comparison.

### Priority B: high-value secondary analyses

1. qPCR-scaled genus-wide differential abundance using MaAsLin 3.
2. qPCR-scale-uncertainty sensitivity using the new ALDEx3 mixed-effects framework.
3. Concordance of relative and qPCR-scaled findings.
4. Replicate repeatability and variance decomposition.
5. Continuous SCFA-response analyses replacing median responder groups.
6. Total SCFA and SCFA-composition analyses.
7. Bray–Curtis versus robust log-ratio community sensitivity.

### Priority C: future-validation analyses

1. Baseline-only prediction with leave-one-donor-out validation.
2. Precision/sample-size simulations for an independent cohort.
3. Direct *Fusicatenibacter*-specific qPCR or ddPCR on retained DNA.
4. Strain-level growth, substrate use, metabolite, co-culture, and isotope-tracing experiments.

### Existing claim gates take precedence

The new qPCR analysis must not displace unfinished primary work:

- The existing butyrate `emmeans` workflow should be generalized to acetate and propionate and exported as one adjusted contrast table.
- Section 5.2 one-way delta ANOVA and pairwise Wilcoxon results do not account for the three repeated carbohydrate conditions within donor. Treat those plots as descriptive or replace their inference with `delta ~ group * carbohydrate + (1 | donor_id)`.
- Donor aggregation fixes well-level pseudo-replication in alpha diversity, but unpaired Wilcoxon tests still discard the repeated donor structure across carbohydrate and time. Fit a donor-aware model for claim-facing alpha results.
- Whole-model PERMANOVA summaries cannot support factor-specific claims. Export term-level R² and p-values and test dispersion before interpretation.
- The Methods mention LEfSe even though the integrated workflow does not implement it. Remove that statement rather than adding another discovery method.
- Reconcile the four absent planned donors and missing 48-hour samples before locking manuscript denominators.

## 4. Phase A: freeze the quantitative data contract

Do not start association modeling until this phase passes.

### 4.1 Build one matched analysis table

Create one row per assay-matched:

```text
donor_id × aliquot_id × carbohydrate × well_repeat × timepoint
```

Required columns:

- `donor_id`;
- `aliquot_id`;
- `well_id`;
- `carbohydrate`;
- `well_repeat`;
- `timepoint`;
- obesity group;
- exact SCFA sample ID;
- exact 16S/qPCR customer label;
- acetate, propionate, butyrate, and total calibrated SCFA;
- ASV read depth;
- *Fusicatenibacter* read count and relative abundance;
- `gene_copies_per_ul`;
- qPCR-scaled *Fusicatenibacter* abundance;
- qPCR Ct and DNA concentration;
- assay-specific missingness and QC flags.

Join SCFA and microbiome on the exact well. Do not average SCFA before this join and then copy the average onto multiple well rows.

### 4.2 Export join diagnostics

Export:

- total rows expected and matched by assay;
- unmatched SCFA IDs;
- unmatched 16S/qPCR IDs;
- duplicate keys;
- label repairs such as repeated periods;
- matched donors, aliquots, wells, and paired 0/48-hour trajectories;
- missingness by group and carbohydrate.

Acceptance:

- no duplicated full join keys;
- all label repairs are explicit and tested;
- every modeled trajectory has both 0 and 48 hours;
- independent donor N is stated in every summary;
- unmatched records are not silently dropped.

### 4.3 qPCR QC

Obtain or document:

- qPCR target and primer sequences;
- amplicon length;
- standard material and standard-curve range;
- efficiency and R²;
- technical replicate rule;
- LOD and LOQ;
- extraction input and elution volume;
- units represented by `gene_copies_per_ul`;
- the derivation and footnote for `genome_copies_per_ul*`;
- handling of Ct values near the assay limit;
- negative and media-control results.

Flag values below LOQ. If LOD/LOQ cannot be recovered, run threshold sensitivities and keep all qPCR-scaled findings explicitly provisional.

## 5. Phase B: focal qPCR-scaled Fusicatenibacter analysis

### 5.1 Primary hypothesis and estimand

Prespecify one focal test:

- exposure: SDC-attributable change in log10 qPCR-scaled *Fusicatenibacter*;
- outcome: SDC-attributable change in butyrate;
- independent unit: donor;
- analysis N: up to 16 donors.

For each well, first calculate 0-to-48-hour change. Average culture-well changes within aliquot, then average A/B aliquots within donor.

Define the substrate-attributable donor-level changes as:

```text
X_d = Δlog10(Fusicatenibacter qPCR-scaled abundance)_SDC
      − Δlog10(Fusicatenibacter qPCR-scaled abundance)_no-carb

Y_d = Δbutyrate_SDC − Δbutyrate_no-carb
```

This difference-in-differences construction controls for time-dependent growth or drift in the no-added-carbohydrate culture.

Primary statistic:

- Spearman correlation between `X_d` and `Y_d`;
- 95% donor-bootstrap confidence interval from 2,000–5,000 replicates;
- exact or donor-level permutation p-value;
- scatter plot with exactly one point per donor and group indicated only for context.

Each bootstrap replicate must resample whole donors with all aliquots, wells, and timepoints intact, then repeat matching, aggregation, qPCR scaling, and model fitting. Prefer percentile/basic intervals because BCa intervals can be unstable with only 16 clusters.

The primary test is prespecified and does not require adjustment for the secondary discovery family. Report the confidence interval even if the p-value is nonsignificant.

### 5.2 Information-preserving sensitivity model

At the aliquot level, retain A/B as replicate observations:

```text
net_delta_butyrate
  ~ net_delta_log10_fusicatenibacter
  + obesity_group
  + (1 | donor_id)
```

Report:

- slope in SCFA units per 10-fold abundance change;
- 95% donor-cluster bootstrap interval;
- marginal and conditional R² where stable;
- residual diagnostics;
- singular-fit status;
- leave-one-donor-out estimates.

Do not include a group-by-*Fusicatenibacter* interaction in the focal model. With only eight donors per obesity group, that interaction is an exploratory sensitivity analysis and should be reported as imprecise.

### 5.3 Zero and detection handling

Lock the zero rule before fitting.

Preferred order:

1. use assay/read-depth information to define a defensible detection floor;
2. analyze log10 abundance with the locked floor;
3. repeat with `log10(1 + abundance)`;
4. repeat donor-level Spearman analysis on untransformed ranks;
5. if zeros remain influential, report detection and positive-abundance components separately.

At 48-hour SDC, only eight of 64 profiles are zero. This is not enough for a reliable condition-specific multivariable prevalence model. Do not force a complex hurdle model solely because MaAsLin 3 can fit one.

### 5.4 Specificity family

Secondary analyses:

1. repeat the net-change analysis for RDC;
2. repeat outcomes for acetate and propionate;
3. analyze total calibrated SCFA;
4. compare relative *Fusicatenibacter* abundance;
5. compare total bacterial qPCR load;
6. estimate the SDC-versus-RDC difference in slopes.

Apply BH correction within this explicitly defined secondary family. Do not infer pathway specificity from “butyrate significant, propionate nonsignificant”; test the relevant slope contrast directly.

### 5.5 Required sensitivity analyses

- mean versus median replicate aggregation;
- raw SDC change versus SDC-minus-no-carb net change;
- gene-copy scaling versus the facility genome-equivalent field;
- all data versus values above LOQ;
- all donors versus leave-one-donor-out;
- healthy-weight and obesity strata shown descriptively;
- exact matched wells versus donor-aggregated table;
- relative versus qPCR-scaled abundance.

An association is considered stable only if its direction does not depend on one donor, one zero rule, or one replicate-aggregation rule.

### 5.6 Interpretation

Allowed:

> Greater SDC-associated change in qPCR-scaled *Fusicatenibacter* abundance was associated with greater SDC-associated butyrate change in this ex vivo dataset.

Not allowed:

- “*Fusicatenibacter* produced butyrate”;
- “*Fusicatenibacter* utilized SDC”;
- “*Fusicatenibacter* mediated the SDC effect”;
- “baseline *Fusicatenibacter* predicts response” unless a separate baseline-only model is validated;
- “absolute abundance measured by *Fusicatenibacter* qPCR.”

## 6. Phase C: total bacterial load

Total load is both a biological outcome and a required check on the focal signal.

Fit:

```text
log10(total 16S gene copies per µL)
  ~ group * carbohydrate * timepoint
  + (1 | donor_id)
  + (1 | aliquot_id)
  + (1 | well_id)
```

If nested variance components are singular, follow a prespecified fallback:

1. donor + aliquot;
2. donor only;
3. donor-level aggregated change model.

Export:

- estimated marginal means;
- direct 0-to-48-hour changes by carbohydrate;
- SDC-versus-RDC and carbohydrate-versus-no-carb change contrasts;
- group modification estimates;
- variance components;
- donor-bootstrap confidence intervals;
- qPCR QC exclusions and sensitivity results.

This determines whether an apparent increase in qPCR-scaled *Fusicatenibacter* is taxon-specific or largely reflects higher total biomass.

## 7. Phase D: qPCR-scaled community-wide analysis

### 7.1 Add MaAsLin 3 as a modern secondary method

MaAsLin 3 is the best candidate new package because it:

- accepts a total-abundance vector from qPCR;
- supports fixed effects, interactions, and random effects;
- separates abundance and prevalence associations;
- supports explicit contrasts;
- can analyze experimental absolute-abundance information without manual TSS inference.

Use the measured total 16S `gene_copies_per_ul` through `unscaled_abundance` with a column named `total`, or pass a precomputed qPCR-scaled table with:

```text
normalization = "NONE"
median_comparison_abundance = FALSE
```

Use a lean model:

```text
~ carbohydrate * timepoint
  + obesity_group
  + (1 | donor_id)
  + (1 | aliquot_id)
```

Do not fit the full group × carbohydrate × time interaction across every genus as the default discovery model. The donor count is too small for stable high-order discovery.

Filter genera by a locked prevalence/abundance rule before modeling. Include read depth in prevalence models. Report abundance and prevalence results separately.

Filter prevalence at the donor level, not the well level. A reasonable starting rule is detection in at least 25% of donors in at least one compared condition, locked before modeling.

### 7.2 Add ALDEx3 as an emerging uncertainty-aware sensitivity

ALDEx3 1.2.0 was released on CRAN on July 15, 2026. It combines Dirichlet-multinomial uncertainty in sequencing counts with uncertainty in the external total-abundance scale and supports mixed-effects models.

For a qPCR-scaled sensitivity analysis:

- supply `s.mu = log2(total 16S gene copies per µL)`;
- derive `s.var` from qPCR technical replication, standard-curve uncertainty, or a prespecified range;
- use the exact `method = "lme4"` engine at this sample size;
- keep filtered taxa in an “other” category so count totals remain coherent;
- retain the same donor/aliquot model and contrast family used by the other methods.

If qPCR variance cannot be estimated, run a range of plausible `s.var` values and show how results change. Do not set scale uncertainty to zero without justification.

Because ALDEx3 is very new and has limited independent benchmarking, use it as a high-value sensitivity analysis rather than the only manuscript method.

### 7.3 Retain ANCOM-BC2 as the relative-abundance comparator

ANCOM-BC2 remains appropriate for compositional repeated-measures differential abundance. Use it to answer the relative/community question, while MaAsLin 3 and ALDEx3 answer quantitative-abundance questions with different treatments of qPCR uncertainty.

For each prespecified contrast, classify genera as:

- concordant direction and supported by both methods;
- qPCR-scaled only;
- relative/compositional only;
- discordant direction;
- too sparse or unstable.

Do not define replication as “significant in one package and nonsignificant in another.” Compare coefficients, confidence intervals, and direction.

### 7.4 Use benchdamic only as a diagnostic

`benchdamic` can summarize distributional fit and method concordance. Use it to justify a small method set, not to run many packages and select whichever gives the smallest q-value.

The manuscript should have:

- one compositional primary method;
- one qPCR-scaled quantitative method;
- one prespecified sensitivity method at most.

## 8. Phase E: other high-value analyses

### 8.1 Replicate repeatability and variance decomposition

Use the corrected hierarchy to quantify:

- donor variance;
- A/B aliquot variance;
- culture-well variance;
- residual/measurement variance;
- intraclass correlations for SCFAs and log10 qPCR load;
- agreement of A/B and R1/R2 or S1/S2 changes.

This directly supports assay reproducibility and determines whether averaging replicates is justified. Report variance components and donor-bootstrap intervals; use Bland–Altman plots only where the compared replicates have the same role.

### 8.2 Continuous response instead of responder groups

Make continuous ΔSCFA the primary interpersonal-response outcome. Retain median responder figures only as descriptive supplement material.

For each focal taxon:

```text
delta_SCFA
  ~ baseline_taxon
  + carbohydrate
  + obesity_group
  + (1 | donor_id)
```

Baseline-only models answer a different, more directionally interpretable question than concurrent 48-hour associations. Keep the covariate count low and perform leave-one-donor-out validation.

### 8.3 Total SCFA and SCFA composition

After facility unit confirmation:

- analyze total calibrated SCFA as a fermentation-output measure;
- analyze acetate, propionate, and butyrate molar fractions as a secondary compositional outcome;
- use log-ratio or multinomial/compositional methods rather than separate unadjusted percentages;
- test whether SDC changes SCFA allocation independently of total production.

Do not combine normalized-only analytes with quantitatively calibrated SCFAs in the total.

### 8.4 Community sensitivity

Complete the already required:

- term-specific PERMANOVA R² and p-values;
- `betadisper`/permutation tests;
- pairwise contrasts with design-matched permutations.

Use separate permutation schemes:

- for carbohydrate or time comparisons, preserve donor pairing and permute condition labels within donor;
- for RDC versus SDC at 48 hours after donor aggregation, use paired within-donor swaps;
- for obesity group, which is constant within donor, aggregate to donor summaries and permute whole donor group labels;
- do not use `strata = donor_id` as the inferential solution for a between-donor obesity effect.

If nuisance covariates or matched-set projection are needed, consider `GUniFrac::adonis3` or LDM as a prespecified sensitivity rather than adding more unconstrained ordinations.

Add one robust log-ratio/Aitchison sensitivity analysis after a locked zero-replacement rule. Bray–Curtis remains the main descriptive analysis; agreement across geometries is more useful than adding more ordinations.

### 8.5 Future-study precision

Use observed donor-level variance and bootstrap distributions to simulate the independent-donor N required to estimate:

- the focal *Fusicatenibacter*–butyrate slope;
- SDC-versus-RDC SCFA contrasts;
- obesity-by-carbohydrate interactions.

Report expected confidence-interval width or assurance, not retrospective “observed power.”

## 9. Reproducible software environment

Do not upgrade the current render environment in place.

The inspected render environment uses:

- R 4.3.3;
- Bioconductor 3.18;
- ANCOMBC 2.4.0;
- Maaslin2 1.18.0;
- no MaAsLin 3;
- no benchdamic;
- `glmmTMB` with a TMB build-version mismatch.

Create a separate locked environment:

- R 4.6.x;
- Bioconductor 3.23;
- `maaslin3` 1.4.x;
- `ALDEx3` 1.2.x;
- `ANCOMBC` 2.14.x;
- `lme4`, `lmerTest`, `nlme`, and `emmeans`;
- `glmmTMB` only if a two-part model is needed, rebuilt against the installed TMB;
- `DHARMa` and/or `performance` for diagnostics;
- `benchdamic` for method-fit and concordance diagnostics;
- `mia` / `TreeSummarizedExperiment` for interoperable data objects;
- `targets` for the analysis DAG;
- `renv` for package locking.

Retain phyloseq as an input compatibility object if useful; a wholesale data-container rewrite is not required to answer the scientific questions.

Record:

- `sessionInfo()`;
- package lockfile;
- random seeds;
- input checksums;
- exact contrast definitions;
- output manifest.

## 10. Proposed file structure

Keep the new analysis outside the 6,000-line integrated Rmd until it passes acceptance.

```text
integrated/additional_analyses/
  _targets.R
  R/
    reconcile_samples_and_missingness.R
    export_primary_scfa_contrasts.R
    alpha_permanova_models.R
    build_matched_scfa_qpcr_table.R
    qpcr_qc.R
    total_load_models.R
    fusicatenibacter_scfa_models.R
    absolute_da_maaslin3.R
    absolute_da_aldex3.R
    concordance_and_sensitivity.R
    export_results.R

integrated/results/additional_analyses_2026-07-19/
  data_audit/
  focal_fusicatenibacter/
  total_load/
  absolute_da/
  diagnostics/
  figures/
  tables/
  manifest.csv
```

After approval, the integrated Rmd should read versioned result tables and figures rather than reimplementing all model logic in report chunks.

## 11. Multiplicity and reporting lock

Prespecify separate families:

1. one focal *Fusicatenibacter*–SDC–butyrate test;
2. focal specificity tests;
3. total-load contrasts;
4. community diversity/PERMANOVA tests;
5. qPCR-scaled genus-wide discovery;
6. relative-abundance genus-wide discovery.

For every claim-facing result, export:

- estimand;
- independent donor N;
- lower-level record counts;
- estimate;
- 95% confidence interval;
- raw p-value;
- adjusted q-value where applicable;
- model formula;
- transformation and zero rule;
- random-effects structure;
- convergence/singularity status;
- donor influence summary.

The focal analysis should not be promoted because it crosses p<0.05. Promotion requires stable direction, acceptable influence diagnostics, correct matching, and biologically restrained interpretation.

## 12. Analyses not recommended for this dataset

Do not add the following to the current manuscript as inferential analyses:

- microbiome co-occurrence networks such as SPIEC-EASI/NetCoMi with N=16 independent donors;
- discovery machine learning or deep learning;
- mediation analysis implying a causal taxon pathway;
- MOFA/DIABLO integration without an independent validation cohort;
- PICRUSt2-based claims of measured metabolic function;
- post hoc subgroup mining across eight donors per obesity group;
- species-level claims from genus-level V3–V4 16S assignments;
- additional responder thresholds selected after inspecting results;
- multiple DA packages used as a significance search.

These methods add complexity without adding reliable independent information at the available donor N.

## 13. Stop/go gates

### Gate 1: data and assay

- [x] Planned-versus-analyzed donors and missing wells are reconciled.
- [x] Adjusted SCFA change contrasts are exported for all three primary analytes.
- [x] Claim-facing alpha models retain donor pairing.
- [x] PERMANOVA term tables and dispersion results are exported.
- [ ] qPCR target, efficiency, standard curve, LOD, and LOQ documented.
- [x] `gene_copies_per_ul` and `genome_copies_per_ul*` definitions reconciled.
- [x] Exact well-level join has no duplicate keys.
- [x] Paired 0/48-hour completeness is exported.
- [ ] Media/negative-control handling is documented.

### Gate 2: focal analysis

- [x] One primary Fusicatenibacter–SDC–butyrate estimand is locked.
- [x] Zero/detection rule is locked before outcome modeling.
- [x] Donor-level N and lower-level records are distinguished.
- [x] Donor bootstrap and leave-one-donor-out results are complete.
- [x] Relative abundance, total load, RDC, and other SCFAs are checked as specificity analyses.

### Gate 3: modern-method expansion

- [ ] Separate `renv` environment is locked.
- [ ] MaAsLin 3 qPCR scaling is verified on a small known example.
- [ ] ALDEx3 scale means and variances are documented and sensitivity-tested.
- [x] Direct contrasts are exported.
- [ ] ANCOM-BC2, MaAsLin 3, and ALDEx3 coefficients are compared, not only q-values.
- [x] Full results, including null and failed fits, are retained.

### Gate 4: manuscript integration

- [ ] Section 5.2 unpaired delta tests are descriptive only or replaced.
- [x] Unimplemented LEfSe language is removed.
- [ ] Existing sample-level correlation and pooled responder p-values are retired.
- [x] Claim audit is updated from versioned outputs.
- [ ] Figures display independent donors honestly.
- [ ] “qPCR-scaled” replaces “direct absolute qPCR” where appropriate.
- [ ] No mechanistic or predictive wording exceeds the design.

## 14. Recommended order of execution

1. Freeze this plan.
2. Reconcile planned/analyzed units and missing samples.
3. Export adjusted carbohydrate-change contrasts for all three primary SCFAs.
4. Replace claim-facing alpha tests and complete PERMANOVA term/dispersion outputs.
5. Obtain qPCR method/QC documentation.
6. Build and validate the exact matched-well table.
7. Run total-load models.
8. Run the focal donor-level *Fusicatenibacter*–butyrate analysis.
9. Run its prespecified sensitivities and specificity family.
10. Quantify replicate repeatability and variance components.
11. Create the locked modern R/Bioconductor environment.
12. Run MaAsLin 3 qPCR-scaled genus-wide analysis.
13. Run the ALDEx3 qPCR-scale-uncertainty sensitivity.
14. Compare with ANCOM-BC2 and existing relative-abundance results.
15. Complete community and SCFA-composition secondary analyses.
16. Refresh the integrated report, claim audit, and manuscript outline.

## 15. Literature and software rationale

Based on articles retrieved from PubMed:

- Quantitative microbiome profiles can differ materially from relative-abundance profiles, but the quantification method itself can introduce substantial technical variability: Galazzo et al., 2020, [doi:10.3389/fcimb.2020.00403](https://doi.org/10.3389/fcimb.2020.00403).
- qPCR-scaled quantitative and relative profiles can produce different abundance patterns: Azarbad et al., 2022, [doi:10.3389/fmicb.2021.798023](https://doi.org/10.3389/fmicb.2021.798023).
- MaAsLin 3 supports prevalence, abundance, repeated measures, and experimental total-abundance information; its methods article was indexed by PubMed as a preprint: Nickols et al., [doi:10.1101/2024.12.13.628459](https://doi.org/10.1101/2024.12.13.628459). The released Bioconductor 3.23 package is version 1.4.0.
- ALDEx3 1.2.0, released July 15, 2026, supports generalized linear/mixed models with sequencing and external-scale uncertainty: [CRAN package](https://cran.r-project.org/package=ALDEx3), [package DOI](https://doi.org/10.32614/CRAN.package.ALDEx3).
- The *F. saccharivorans* type-strain paper reported acetate, lactate, formate, and succinate as glucose fermentation products and does not establish direct butyrate production: Takada et al., 2013, [doi:10.1099/ijs.0.045823-0](https://doi.org/10.1099/ijs.0.045823-0).

ANCOM-BC2 remains the established repeated-measures compositional comparator: Lin and Peddada, 2024, [doi:10.1038/s41592-023-02092-7](https://doi.org/10.1038/s41592-023-02092-7).
