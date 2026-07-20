# Gut Microbes manuscript claim audit

**Project:** Abbott–Lurie ZC07, Project 2  
**Date:** July 19, 2026  
**Updated:** After incorporation of Sections 5.3.1 and 5.4.2 from the donor-aware integrated re-knit

**Primary evidence source:** `render/integrated-scfa-microbiome.html`

**Purpose:** Lock manuscript claims to the current donor-aware analyses and identify remaining analytical, methods, and reporting gates.

## Executive decision

The current dataset supports an ex vivo functional-phenotyping paper centered on analyte-specific SCFA trajectories and repeatable donor heterogeneity. The new donor-aware 48-hour endpoint analysis additionally supports an omnibus carbohydrate-condition association with absolute acetate, propionate, butyrate, and succinate concentrations after BH adjustment. It does not support SDC superiority over RDC for acetate, propionate, or butyrate; equivalence of obesity groups; a global carbohydrate effect on 48-hour community structure; a significant donor-level *Fusicatenibacter*–butyrate association; a stable responder phenotype; prediction; or microbial mechanism.

The strongest defensible conclusion is:

> Ex vivo SCFA accumulation by adolescent fecal microbiota varied with incubation time, carbohydrate condition, and donor identity. At 48 hours, carbohydrate condition was associated with endpoint acetate, propionate, butyrate, and succinate concentrations in donor-aware models. This omnibus endpoint result does not identify a superior formulation: adjusted changes under SDC did not differ significantly from RDC for the three primary SCFAs. Obesity status did not significantly modify SCFA changes, although equivalence was not established. Time and obesity group were associated with whole-community structure, while integrated taxon–SCFA findings were limited and exploratory.

The manuscript can proceed to a claim-locked Results shell. Numerical units, exclusions, cohort details, fermentation methods, and several render defects must be resolved before the Abstract and final Results prose are frozen.

## Evidence precedence

The current HTML contains both new claim-lock sections and older summary prose. Use evidence in this order:

1. The donor-aware 48-hour endpoint term summary in Section 5.3.1 for omnibus group, carbohydrate, and group×carbohydrate inference.
2. The repaired live nested-LMM contrast in Section 5.4.2 for the prespecified butyrate SDC-versus-RDC change estimand.
3. Versioned additional-analysis exports loaded in Sections 6.7.1–6.7.4, except where a newer live result explicitly supersedes an older export row.
4. Other donor-aware SCFA results in Sections 5.2, 5.4, and 5.5.
5. Term-level community results and dispersion tests in Section 6.7.2.
6. Current MaAsLin3, ALDEx3, and ANCOM-BC2 outputs in Sections 6.5.6–6.5.14.
7. Descriptive sample-level and median-split panels only when explicitly labeled exploratory.

Do not use Section 8 as an independent result source where it conflicts with earlier claim locks. It retains superseded alpha-, beta-, and *Fusicatenibacter* statements.

## Evidence classes

- **Supported:** Directly tested with the current donor-aware model and interpretable from the reported estimate and uncertainty.
- **Supported with qualification:** The reported result is valid, but its wording must reflect effect size, multiplicity, model scope, or design limitations.
- **Exploratory:** A reproducible analysis produced the result, but the threshold, multiplicity, design, or lack of independent validation precludes a primary claim.
- **Unsupported:** The current result contradicts the claim or the required test was not performed.
- **Retired:** A prior or descriptive result has been explicitly replaced by a more appropriate claim-facing analysis.

## Claim-by-claim audit

| Proposed claim | Classification | Current evidence | Manuscript-safe wording |
|---|---|---|---|
| Incubation time and carbohydrate condition shape SCFA trajectories. | Supported | Acetate, propionate, and butyrate increased from 0 to 48 hours in all three conditions. Mixed models reported carbohydrate and selected carbohydrate-by-time effects, while direct contrasts showed analyte-specific patterns. | Incubation time and carbohydrate condition were associated with analyte-specific ex vivo SCFA trajectories. |
| Carbohydrate condition is associated with absolute SCFA concentrations at 48 hours. | Supported with scope qualification | In donor-averaged 48-hour models of `concentration ~ group × carbohydrate + (1|donor)`, the omnibus carbohydrate term was BH-significant for acetate, propionate, butyrate, and succinate. It was not significant for 5-aminovalerate; no group or group×carbohydrate term was BH-significant. | At 48 hours, carbohydrate condition was associated with endpoint acetate, propionate, butyrate, and succinate concentrations in donor-aware models. |
| Both RDC and SDC increase acetate, propionate, and butyrate relative to no added carbohydrate. | Unsupported as a general inferential claim | All conditions, including no added carbohydrate, accumulated SCFAs. Most adjusted active-carbohydrate-versus-control differences in change were nonsignificant; RDC had a smaller propionate change than control. | SCFAs increased under RDC and SDC, but adjusted comparisons with no added carbohydrate were not uniformly positive or significant. |
| SDC produces more acetate, propionate, or butyrate than RDC. | Unsupported | SDC-minus-RDC differences in change were acetate 2.89 (95% CI −2.06 to 7.84; p=0.540) and propionate 1.06 (−0.32 to 2.44; p=0.230) in the earlier pairwise export. The repaired prespecified butyrate contrast gave mean changes of 7.95 μM under SDC and 6.92 μM under RDC and an adjusted difference of 1.03 μM (−0.85 to 2.91; p=0.396). | Adjusted 0-to-48-hour changes did not differ significantly between SDC and RDC for any primary SCFA. |
| SDC increases butyrate relative to no added carbohydrate. | Provisional; source conflict | The older export gives a difference of 2.62 (95% CI −0.001 to 5.23; Tukey-adjusted p=0.0502), while Section 5.6 calls the contrast significant. The repaired Section 5.4.2 text supplies only the SDC-versus-RDC row. | Do not claim an SDC-versus-control difference until all three live butyrate pairwise rows are exported and reconciled. |
| Obesity impairs or preserves fermentation capacity. | Unsupported | None of six obesity-minus-healthy-weight differences in 0-to-48-hour change was significant. No externally justified equivalence margin was available. | SCFA changes did not differ significantly by obesity status in this cohort; equivalence was not evaluated. |
| Donors differ in generalized SCFA responsiveness. | Supported | Donor ICCs were 0.63–0.84; RDC–SDC rank correlations were 0.83–0.96 across the three primary analytes. | Donors showed repeatable high-to-low SCFA response ordering across carbohydrate conditions. |
| Donors have stable substrate-specific response phenotypes. | Unsupported | Strong cross-condition correlations indicate generalized responsiveness. A donor-by-carbohydrate variance component and independent replication were not reported. | The data support generalized donor response heterogeneity, not stable substrate-specific phenotypes. |
| Obesity-group SCFA response dispersion differs. | No detected difference; imprecise | Obesity-to-healthy-weight SD ratios ranged from 0.66 to 1.23; no Brown–Forsythe comparison survived BH adjustment (minimum q=0.995), with 8 donors per group. | No obesity-group difference in SCFA response dispersion was detected; precision was limited. |
| Healthy-weight donors have higher alpha diversity at 48-hour no-carb. | Supported with multiplicity qualification | In donor-aware mixed-model contrasts, obesity-minus-healthy-weight estimates were −0.504 for Shannon (p=0.00425) and −42.8 for Observed richness (p=0.00241); Simpson was nonsignificant. A post hoc BH calculation across all 18 exported group contrasts gives q≈0.038 for Shannon and Observed, but adjusted q-values are not yet versioned in the export. | At 48 hours without added carbohydrate, Shannon diversity and observed richness were lower in the obesity group in donor-aware models. State the locked multiplicity family and adjusted q-values. |
| Carbohydrate, time, group, and group-by-carbohydrate each affect beta diversity. | Unsupported | In the donor-aggregated all-timepoint model, time R²=0.2390, p=0.0001 and group R²=0.0264, p=0.0037; carbohydrate p=0.665 and group×carbohydrate p=1.000. At 48 hours, group R²=0.0555, p=0.0027; carbohydrate p=0.4085 and interaction p=1.000. | Culture time and obesity group were associated with Bray–Curtis structure; carbohydrate and group-by-carbohydrate terms were not significant in the donor-aggregated models. |
| The group PERMANOVA is solely a dispersion artifact. | Not supported by the current dispersion test | At 48 hours, betadisper p=0.0942 for group and p=0.8759 for carbohydrate. | No significant difference in multivariate dispersion was detected, although the group result was imprecise. |
| Multiple genera robustly distinguish RDC from SDC at 48 hours. | Unsupported | Focused ANCOM-BC2 reported one genus at exploratory q<0.1. The exact head-to-head coefficient, interval, and q-value are not yet in a versioned manuscript table. | One relative-abundance genus met the exploratory q<0.1 threshold in the focused 48-hour comparison. |
| Relative *Fusicatenibacter* abundance is higher under SDC than RDC at 48 hours. | Exploratory | Focused ANCOM-BC2 identified *Fusicatenibacter* as the single q<0.1 genus, higher under SDC. The result is relative, uses an exploratory threshold, and lacks a versioned exact contrast export. | Relative *Fusicatenibacter* abundance was higher under SDC than RDC at 48 hours in an exploratory ANCOM-BC2 analysis. |
| SDC expands qPCR-scaled *Fusicatenibacter*. | Unsupported | MaAsLin3 SDC×48-hour coefficient 3.78, q=0.125; ALDEx3 coefficient 4.31, adjusted p=0.581. The relative ANCOM reference-coded SDC_48H coefficient was not supportive of an absolute expansion claim. | qPCR-scaled analyses did not establish a significant SDC-associated expansion of *Fusicatenibacter*. |
| Differential-abundance methods are broadly concordant. | Unsupported as a broad claim | Among 182 focus rows with all three coefficients, 85 (47%) were direction-concordant. For SDC×48 hours, 16 of 45 were three-way concordant. | Relative, qPCR-scaled, and scale-uncertainty methods provided partially discordant views and should be reported as complementary estimands. |
| *Bifidobacterium* increases with both carbohydrates, more under SDC. | Unsupported in this form | Current MaAsLin3 found an exploratory SDC×48-hour coefficient of 3.28 (q=0.0745); RDC×48-hour q=0.129. ANCOM-BC2 relative coefficients were nonsignificant, and no direct SDC-versus-RDC qPCR-scaled contrast was reported. | qPCR-scaled *Bifidobacterium* showed an exploratory SDC×48-hour interaction in MaAsLin3; the analysis does not show significant increases under both carbohydrates or SDC superiority. |
| Relative *Fusicatenibacter* correlates significantly with butyrate under SDC. | Retired | The descriptive well-level 48-hour SDC analysis gave R²=0.0718, nominal p=0.0459, Spearman ρ=0.365, N=56, but treats nested wells as independent and is explicitly retired in Section 7.4. | Do not report the sample-level p-value as claim-facing evidence. |
| Donor-level qPCR-scaled *Fusicatenibacter* change associates with SDC-attributable butyrate change. | Exploratory, nonsignificant | Exact well-matched donor analysis: Spearman ρ=0.477, permutation p=0.088, bootstrap 95% CI −0.239 to 0.892, N=14. All specificity-family BH q-values were ≥0.949. | The donor-level association was positive but imprecise and did not reach statistical significance. |
| The *Fusicatenibacter* signal is butyrate- or pathway-specific. | Unsupported | Specificity-family tests for RDC butyrate, SDC relative abundance, total load, acetate, propionate, and total SCFA were all nonsignificant after BH adjustment. Significant-versus-nonsignificant comparisons do not test effect differences. | No pathway-specific association was established. |
| qPCR-scaled *Fusicatenibacter* is higher in butyrate responders. | Retired as a claim-facing test; descriptive only | Pooled well-level Wilcoxon p=0.00123, with 63 records per response group. It dichotomizes the outcome, pools RDC and SDC, and retains multiple observations per donor. | Pooled descriptive abundances differed across median-defined strata, but the donor-level continuous analysis is the inferential result. |
| *Bifidobacterium* separates SCFA response strata. | Unsupported | Pooled qPCR-scaled Wilcoxon p=0.103 for butyrate strata and p=0.531 for propionate strata. | qPCR-scaled *Bifidobacterium* did not differ significantly by median-defined response stratum. |
| Butyrate response strata have distinct community structure. | Unsupported by the current PERMANOVA | The donor-stratified responder×carbohydrate whole model gave R²=0.0875, p=0.108. Pooled alpha comparisons were nominally significant but remain concurrent and median-defined. | Median-defined butyrate strata showed exploratory alpha-diversity differences but no significant whole-community PERMANOVA result. |
| Propionate response strata have distinct community structure. | Unsupported | Alpha comparisons were nonsignificant after adjustment; responder×carbohydrate PERMANOVA R²=0.0500, p=0.192. | No alpha- or beta-diversity difference was detected for median-defined propionate strata. |
| Butyrate and propionate responder pathways are independent. | Unsupported | The symmetric overlap table contains 12 both higher, 12 both lower, and 12 in each discordant cell across 48 donor×carbohydrate records. Median splitting forces balanced groups, and repeated donor records are not independent. | The response classifications showed incomplete descriptive overlap; pathway independence was not tested. |
| Median-defined responders constitute a stable phenotype. | Unsupported | A/B aliquots fell on the same side of the donor-derived threshold for only 7–12 of 16 donors depending on analyte and carbohydrate. | Median-defined higher- and lower-response strata are sample-dependent exploratory groupings. |
| *Fusicatenibacter* uses SDC or produces the observed butyrate. | Unsupported mechanistically | No isolate growth, substrate disappearance, isotope flux, metatranscriptomic, or co-culture data were collected. | The genus is a candidate for mechanistic follow-up; the current data do not establish SDC utilization or butyrate production. |
| The ex vivo assay or microbiome predicts dietary response. | Unsupported | Microbiome and SCFA features are mostly concurrent, the outcome is defined within the same experiment, and no baseline-only held-out validation or in vivo outcome was analyzed. | The study generates hypotheses for future response prediction; it does not validate a predictor. |

## Current numerical anchors

### Donor-aware SCFA concentrations at the 48-hour endpoint

Section 5.3.1 averages A/B aliquots and culture wells to one value per donor × carbohydrate × analyte, then fits `concentration_48H ~ group × carbohydrate + (1|donor)` separately by analyte. Type III term p-values are BH-adjusted across the 15 analyte × model-term tests.

| Analyte | Carbohydrate F(df) | Raw p | BH q | Claim status |
|---|---:|---:|---:|---|
| Propionate | 25.1612 (2, 28) | 5.57×10⁻⁷ | 1.67×10⁻⁶ | Supported omnibus endpoint association |
| Acetate | 16.1075 (2, 28) | 2.21×10⁻⁵ | 6.63×10⁻⁵ | Supported omnibus endpoint association |
| Butyrate | 8.5831 (2, 28) | 0.00124 | 0.00371 | Supported omnibus endpoint association |
| Succinate | 5.4613 (2, 28) | 0.00994 | 0.0298 | Supported omnibus endpoint association |
| 5-aminovalerate | 1.4082 (2, 28) | 0.261 | 0.784 | No BH-significant endpoint term |

No obesity-group or group×carbohydrate term was BH-significant. The closest nonsignificant interaction was acetate, F(2,28)=2.9422, raw p=0.0692, q=0.104. The complete 15-row table is embedded in the rendered HTML and exported as `integrated/results/additional_analyses_2026-07-19/scfa_48h_endpoint_type3_anova.csv`.

This endpoint result answers whether final 48-hour concentrations differ somewhere across carbohydrate conditions. It does **not** identify which pair differs, show that either active carbohydrate exceeds no added carbohydrate, estimate a 0-to-48-hour response, or establish SDC superiority over RDC. Those claims require the adjusted pairwise change contrasts below. Descriptive Wilcoxon stars in Section 5.3.1 are not claim-facing evidence.

### Primary carbohydrate differences in 0-to-48-hour change

Values use the concentration unit labeled µM in the workflow, pending facility confirmation. The repaired Section 5.4.2 live nested-LMM result is authoritative for the prespecified butyrate SDC-versus-RDC contrast and supersedes the older row in `carbohydrate_change_contrasts.csv`. Other rows below remain from the earlier exported family and should be refreshed from the same live model before the manuscript table is locked.

| Analyte | Contrast | Estimate | 95% CI | Tukey-adjusted p |
|---|---|---:|---:|---:|
| Acetate | RDC − no added carbohydrate | 3.00 | −3.13 to 9.12 | 0.726 |
| Acetate | SDC − no added carbohydrate | 5.89 | −0.23 to 12.01 | 0.066 |
| Acetate | SDC − RDC | 2.89 | −2.06 to 7.84 | 0.540 |
| Propionate | RDC − no added carbohydrate | −1.86 | −3.57 to −0.15 | 0.025 |
| Propionate | SDC − no added carbohydrate | −0.80 | −2.51 to 0.91 | 0.764 |
| Propionate | SDC − RDC | 1.06 | −0.32 to 2.44 | 0.230 |
| Butyrate | RDC − no added carbohydrate | 1.58 | −1.03 to 4.20 | 0.500 |
| Butyrate | SDC − no added carbohydrate | 2.62 | −0.001 to 5.23 | 0.0502 |
| Butyrate | SDC − RDC | 1.03 | −0.85 to 2.91 | 0.396 |

The repaired butyrate model estimated mean 0-to-48-hour changes of 7.95 μM under SDC and 6.92 μM under RDC. The 1.03-μM adjusted difference was nonsignificant, so the analysis does not demonstrate a greater butyrate response under SDC.

In the earlier general export, the only active-carbohydrate-versus-control difference meeting p<0.05 was a smaller propionate increase under RDC than under no added carbohydrate. Background SCFA accumulation in the control medium argues against a simple uniformly increased-output narrative. Because the repaired live butyrate row differs from that export, refresh all three live butyrate pairwise rows before locking the SDC-versus-control claim.

### Obesity-group differences in change

The prespecified estimand is:

`(obesity 48 h − obesity 0 h) − (healthy-weight 48 h − healthy-weight 0 h)`

| Analyte | Carbohydrate | Estimate | 95% CI |
|---|---|---:|---:|
| Acetate | RDC | −0.09 | −5.17 to 5.00 |
| Acetate | SDC | 3.43 | −1.66 to 8.51 |
| Propionate | RDC | 0.74 | −0.70 to 2.18 |
| Propionate | SDC | 0.38 | −1.06 to 1.82 |
| Butyrate | RDC | −0.76 | −2.96 to 1.44 |
| Butyrate | SDC | −1.19 | −3.39 to 1.02 |

None differed significantly from zero. Donor-level bootstrap intervals gave the same qualitative result. Residual diagnostics found non-normality for all three analytes and heteroscedasticity for propionate and butyrate; retain bootstrap sensitivity intervals in the supplement.

No external biological threshold or project-specific assay-precision bound was available. The equivalence-margin file remains unlocked, so TOST p-values and equivalence verdicts should not be reported. Post hoc precision frontiers are not equivalence margins.

### Donor heterogeneity and repeatability

| Analyte | Donor ICC (bootstrap 95% CI) | RDC–SDC Spearman ρ (bootstrap 95% CI) |
|---|---:|---:|
| Acetate | 0.84 (0.49 to 0.95) | 0.96 (0.83 to 1.00) |
| Propionate | 0.84 (0.66 to 0.90) | 0.85 (0.51 to 0.96) |
| Butyrate | 0.63 (0.41 to 0.76) | 0.83 (0.51 to 0.97) |

These statistics now support a quantitative interpersonal-variation claim. They point to generalized high-to-low responsiveness rather than donor-specific preference for SDC or RDC.

### Alpha diversity

The current claim-facing export uses donor-aware mixed models. At 48 hours under no added carbohydrate, obesity-minus-healthy-weight contrasts were:

- Shannon: −0.504 (95% CI −0.843 to −0.165; p=0.00425);
- Observed richness: −42.8 (−69.8 to −15.8; p=0.00241);
- Simpson: −0.0349 (−0.0929 to 0.0232; p=0.234).

A post hoc BH adjustment across all 18 exported group contrasts yields q≈0.0383 for Shannon and Observed richness. Because the export currently contains only raw p-values, add the prespecified adjustment and q-value columns to the versioned table before manuscript citation. Do not use diversity as a health or improvement endpoint.

### Bray–Curtis community structure

| Model | Term | R² | p |
|---|---|---:|---:|
| All timepoints, donor aggregated | Group | 0.0264 | 0.0037 |
| All timepoints, donor aggregated | Carbohydrate | 0.0140 | 0.6651 |
| All timepoints, donor aggregated | Time | 0.2390 | 0.0001 |
| All timepoints, donor aggregated | Group × carbohydrate | 0.0041 | 1.0000 |
| 48 hours, donor aggregated | Group | 0.0555 | 0.0027 |
| 48 hours, donor aggregated | Carbohydrate | 0.0435 | 0.4085 |
| 48 hours, donor aggregated | Group × carbohydrate | 0.0095 | 1.0000 |

At 48 hours, betadisper p=0.8759 for carbohydrate and p=0.0942 for group. Whole-model PERMANOVA values from older sections must not be used to assign significance to every term.

### Focal *Fusicatenibacter* results

- Focused relative ANCOM-BC2 at 48 hours: one q<0.1 genus, *Fusicatenibacter*, higher under SDC than RDC. Exact coefficient, interval, and q-value still need a versioned export.
- MaAsLin3 qPCR-scaled SDC×48-hour coefficient: 3.777, q=0.125; not significant.
- ALDEx3 SDC×48-hour coefficient: 4.308, adjusted p=0.581; not significant.
- Donor-level SDC-attributable qPCR-scaled change versus butyrate change: ρ=0.477, permutation p=0.088, bootstrap 95% CI −0.239 to 0.892, N=14; not significant.
- All prespecified specificity-family BH q-values were ≥0.949.

## Experimental-unit and sample audit

| Level | Current evidence | Interpretation |
|---|---|---|
| Planned cohort | 20 numeric participants / 40 A/B aliquot labels in facility metadata | Protocol or submission set, not analyzed SCFA N |
| SCFA donors | 16 numeric donors: 8 healthy-weight and 8 obesity | Independent biological unit for group inference |
| SCFA aliquots | 32 A/B stool splits | Nested within donor; not independent people |
| SCFA culture wells | 160 in the corrected nested workflow | Nested within aliquot × carbohydrate and linked across 0/48 hours |
| Microbiome profiles | 292 phyloseq samples, all matched to qPCR | Sample-level records; donor-aware aggregation or clustering required |
| ANCOM-BC2 analysis | 290 samples after analysis-specific cleaning; 20 genera | Report analysis-specific exclusions |
| 48-hour focused ANCOM | 126 samples: 62 RDC and 64 SDC | Nested records, not 126 participants |
| Responder overlap | 48 donor × carbohydrate records | Repeated records from 16 donors, not 48 people |

The manuscript must explain why four planned numeric donors—eight A/B labels—are absent from the SCFA export and why two obesity no-added-carbohydrate 48-hour observations are missing.

## Methods and terminology audit

### Timepoints

The analyzed files and current integrated report contain 0- and 48-hour observations. Do not include 72 hours unless a separate dataset is supplied and analyzed.

### SCFA assay

Locally documented elements include PFBBr derivatization, negative-chemical-ionization GC-MS, and authentic-standard calibration for acetate, propionate, butyrate, and succinate. Resolve before drafting:

- exported concentration unit and dilution mapping;
- analyte-specific CV, LOD, and LOQ;
- Agilent 8890 overview versus 7890A/5975C detailed method;
- 500 mg submitted culture material versus 100 mg/mL facility extraction description; and
- quantitatively calibrated concentrations versus normalized-abundance analytes.

Use “SCFA accumulation” or “net output.” Concentration change does not measure production flux.

### Fermentation design

The repository still lacks manuscript-ready documentation for exact RDC/SDC composition, source, lot, dose, medium, inoculum, anaerobic atmosphere, incubation temperature, vessel volume, mixing, pH, plate allocation, batch handling, and background carbohydrate in the no-added-carbohydrate condition.

The design should be described as direct ex vivo exposure to formulations classified by expected host digestibility. The experiment itself does not measure digestion.

### 16S and qPCR

The current workflow uses Zymo Bac16S V3–V4 data, DADA2-derived ASVs, genus aggregation, phyloseq, ANCOM-BC2, MaAsLin3, and ALDEx3. Add primer sequences, filtering parameters, taxonomy database/version, contaminant handling, depth exclusions, accession, qPCR target, standard curve, efficiency, and LOD/LOQ.

Define the focal quantitative measure as genus relative abundance × total bacterial 16S target load per µL. Use “qPCR-scaled genus abundance” or “estimated 16S target equivalents,” not genus-specific qPCR, cells, or unqualified absolute abundance.

### Terminology

- Confirm that FDC in project descriptions maps to RDC in the executed analysis.
- Confirm whether GFA was omitted, renamed, or represented by an analyzed condition.
- Standardize Control/Case to healthy-weight/obesity after verifying the clinical definition.
- Use “higher-response” and “lower-response strata,” not validated responders.

## Render discrepancy ledger

The current HTML is a useful source only when its internal precedence is respected.

1. **Section 5.4.2 is repaired and now authoritative.** It reports mean changes of 7.95 μM under SDC and 6.92 μM under RDC and an adjusted difference of 1.03 μM (95% CI −0.85 to 2.91; p=0.396).
2. **Section 5.3.1 adds a distinct, supported endpoint claim.** Its BH-significant omnibus carbohydrate terms apply to 48-hour absolute concentrations of acetate, propionate, butyrate, and succinate. They do not resolve pairwise carbohydrate differences or changes from baseline.
3. **The Section 5.3.1 model table is now static and versioned.** The HTML embeds all 15 Type III rows, and the exact F statistics, degrees of freedom, raw p-values, and BH q-values are exported to `scfa_48h_endpoint_type3_anova.csv`.
4. **The older pairwise CSV is stale for the prespecified butyrate contrast.** It reports the same point estimate but a different interval and p-value. Refresh or version the export so downstream tables use the repaired result.
5. **The SDC-versus-control butyrate claim remains inconsistent.** The older export gives p=0.0502, while Section 5.6 calls the contrast significant. Export all three live Section 5.4.2 pairwise rows before locking this secondary claim.
6. **Sections 7.2 and 7.3 contain parse errors.** No continuous ANCOM-BC2 butyrate or propionate discovery result should be claimed from these sections.
7. **Section 8.2 repeats superseded beta-diversity claims.** Carbohydrate and group×carbohydrate were not significant in the donor-aggregated term tables.
8. **Section 8.2 does not reflect the donor-aware alpha export cleanly.** Use the mixed-model estimates and add a locked multiplicity adjustment.
9. **Section 8.3 repeats the retired sample-level *Fusicatenibacter* correlation.** The nominal p=0.0459 is not claim-facing; the donor-level permutation p=0.088 and bootstrap interval are authoritative.
10. **The focused relative ANCOM result is not versioned as a manuscript table.** Export the exact coefficient, interval, and q-value.
11. **Several exploratory pooled responder plots remain visible after claim retirement.** Their presence does not restore inferential status.

## Completed analysis gates

- [x] Correct numeric donor → A/B aliquot → culture-well hierarchy implemented.
- [x] Donor-aware 48-hour SCFA endpoint models run with BH adjustment across analyte × model terms.
- [x] Complete Section 5.3.1 Type III endpoint table embedded in HTML and exported as a versioned CSV.
- [x] Repaired live nested-LMM SDC-versus-RDC butyrate contrast rendered with adjusted means, interval, and p-value.
- [x] Earlier adjusted RDC-versus-SDC and active-carbohydrate-versus-control SCFA contrast family exported.
- [x] Six obesity-group difference-in-change contrasts and bootstrap sensitivities completed.
- [x] Donor ICCs and RDC–SDC rank correlations reported.
- [x] Donor-aware alpha mixed-model contrasts exported.
- [x] Term-level donor-aggregated PERMANOVA exported.
- [x] Betadisper results exported.
- [x] MaAsLin3 qPCR-scaled and ALDEx3 scale-uncertainty screens completed.
- [x] Exact well-matched donor-level *Fusicatenibacter*–butyrate analysis completed.
- [x] Sample-level correlation and pooled responder Wilcoxon results explicitly retired from claim-facing use.

## Remaining gates

### Blocking final Results and Abstract

1. Confirm SCFA concentration units, dilution mapping, CV, LOD, and LOQ.
2. Reconcile the four absent planned donors and missing SCFA observations.
3. Supply cohort, eligibility, obesity-definition, demographic, and ethics details.
4. Supply the complete fermentation and substrate protocol.
5. Refresh the versioned SCFA pairwise table from the repaired live model and reconcile the remaining butyrate contrasts.
6. Repair the failed continuous ANCOM chunks, then re-knit.
7. Export the focused RDC-versus-SDC ANCOM-BC2 coefficient, interval, and q-value.
8. Add the locked alpha-diversity multiplicity family and q-values to the versioned export.

### Required only for stronger optional claims

1. Add a donor-by-carbohydrate variance component before claiming substrate-specific donor phenotypes.
2. Use baseline-only predictors and held-out validation before prediction language.
3. Fit direct response-type interactions before claiming butyrate specificity.
4. Perform isolate, substrate-use, flux, or co-culture experiments before mechanistic attribution.

Equivalence margins are not a drafting gate if the manuscript states that equivalence was not evaluated. They become necessary only if the study team wants an equivalence or non-inferiority claim.

## Language guardrails

Use:

- “associated with”
- “net accumulation”
- “did not differ significantly”
- “equivalence was not evaluated”
- “repeatable donor response ordering”
- “qPCR-scaled genus abundance”
- “relative *Fusicatenibacter* signal”
- “positive but nonsignificant”
- “median-defined response strata”
- “hypothesis-generating”

Avoid:

- “SDC is superior”
- “both carbohydrates increased every SCFA versus control”
- “obesity did not affect”
- “preserved,” “equivalent,” or “unaffected capacity”
- “carbohydrate drove community structure at 48 hours”
- “absolute *Fusicatenibacter* expansion”
- “significant *Fusicatenibacter*–butyrate correlation”
- “pathway-specific”
- “SDC-utilizing butyrate producer”
- “responder phenotype”
- “predicts response”
- “companion diagnostic”

## Recommended title after the updated audit

**Preferred**

> Incubation Time, Carbohydrate Condition, and Donor Identity Shape Ex Vivo Short-Chain Fatty Acid Accumulation by Adolescent Fecal Microbiota

**Alternative**

> Interpersonal Variation in Ex Vivo Short-Chain Fatty Acid Responses of Adolescent Fecal Microbiota to Carbohydrate Challenge

Do not put SDC superiority, obesity independence, responder phenotypes, prediction, or *Fusicatenibacter* in the title.

## Source ledger

- `render/integrated-scfa-microbiome.html`: current integrated render; use claim-lock sections rather than conflicting final summary prose.
- `render/integrated-scfa-microbiome.html`, Section 5.3.1: donor-aware 48-hour endpoint Type III table and omnibus term summary.
- `render/integrated-scfa-microbiome.html`, Section 5.4.2: authoritative repaired live nested-LMM butyrate SDC-versus-RDC contrast.
- `integrated/results/additional_analyses_2026-07-19/scfa_48h_endpoint_type3_anova.csv`: exact 15-row endpoint term export with raw p-values and BH q-values.
- `integrated/results/additional_analyses_2026-07-19/scfa_contrasts/carbohydrate_change_contrasts.csv`: earlier general SCFA change contrasts; its butyrate SDC-versus-RDC row is superseded and must be refreshed.
- `integrated/results/additional_analyses_2026-07-19/alpha_permanova/`: donor-aware alpha, term-level PERMANOVA, and betadisper exports.
- `integrated/results/additional_analyses_2026-07-19/focal_fusicatenibacter/`: donor-level focal association, specificity family, bootstrap, and influence outputs.
- `integrated/results/additional_analyses_2026-07-19/absolute_da/`: MaAsLin3, ALDEx3, ANCOM-BC2 reference-coded results, and concordance tables.
- `integrated/results/additional_analyses_2026-07-19/data_audit/`: donor, trajectory, missingness, and join audit.
- `scfa_metabolomics/obesity_equivalence_analysis.py` and `scfa_metabolomics/results/obesity_group_scfa_*.csv`: obesity-group contrasts, diagnostics, and bootstrap sensitivities.
- `scfa_metabolomics/equivalence_margins.csv`: unlocked margin specification; no equivalence verdict.
- `docs/scfa-obesity-equivalence-margin-spec.2026-07-19.md`: equivalence rationale and precision-frontier interpretation.
- `manuscript/manuscript-outline.md`: manuscript structure synchronized to this updated audit.

## Current go/no-go decision

**Go:** claim-locked manuscript structure, Introduction, provisional Methods, SCFA Results shell, community Results shell, figures, and literature mapping.

**Conditional go:** final numerical Results after unit confirmation, exclusion reconciliation, alpha q-value export, focused ANCOM export, and clean re-knit.

**No-go:** SDC superiority, equivalence, mechanistic *Fusicatenibacter*, validated responder, predictive, clinical-benefit, or personalized-treatment claims.
