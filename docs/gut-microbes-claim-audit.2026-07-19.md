# Gut Microbes manuscript claim audit

**Project:** Abbott–Lurie ZC07, Project 2  
**Date:** July 19, 2026  
**Purpose:** Determine which proposed manuscript claims are supported by the current analysis, which require qualification, and which require reanalysis or additional source information before drafting a *Gut Microbes* submission.

## Executive decision

The current dataset supports an ex vivo phenotyping paper centered on incubation time, carbohydrate condition, limited modification of SCFA trajectories by obesity status, and marked interpersonal variation. It does not yet support the proposed mechanistic title or a validated responder-prediction claim.

The manuscript can proceed to a structured shell, literature mapping, and provisional Methods. Claim-locked Results and Abstract drafting should wait until the donor/replicate structure is resolved and the primary RDC-versus-SDC contrasts are rerun.

The safest current conclusion is:

> In this ex vivo system, incubation time and carbohydrate condition were associated with SCFA production, while obesity status did not significantly modify SCFA trajectories. Community and taxon associations, including the 48-hour absolute-abundance signal for *Fusicatenibacter*, are exploratory and require condition-specific and independent validation.

## Evidence classes

- **Supported:** Directly tested by an appropriate reported model and interpretable from the current output.
- **Supported with qualification:** The direction is plausible or partly demonstrated, but wording must reflect missing contrasts, multiplicity, uncertainty, or design constraints.
- **Exploratory:** A real analysis produced the result, but the design does not establish prediction, mechanism, or causality.
- **Unsupported:** The current output does not provide the required test, or the reported result contradicts the claim.

## Claim-by-claim audit

| Proposed claim | Classification | Evidence | Manuscript-safe wording |
|---|---|---|---|
| Carbohydrate condition and time are the primary drivers of SCFA production. | Supported | Mixed models found strong time effects for all five analytes and carbohydrate effects for acetate, butyrate, and propionate. BH-adjusted time p-values were approximately 4.9×10⁻⁴⁶ for acetate, 5.4×10⁻³² for propionate, 3.1×10⁻²⁷ for butyrate, 2.2×10⁻²⁴ for 5-aminovalerate, and 1.7×10⁻⁸ for succinate. Carbohydrate main-effect adjusted p-values were 0.011 for acetate, 0.042 for butyrate, and 0.0084 for propionate. | Incubation time and carbohydrate condition were associated with ex vivo SCFA concentrations. |
| Both RDC and SDC increase acetate, butyrate, and propionate relative to no added carbohydrate. | Supported with qualification | The global carbohydrate effect on 48-hour deltas was significant for acetate (p=0.0479), butyrate (p=0.0420), and propionate (p=0.0209). The rendered report does not provide adjusted pairwise no-carb versus RDC and no-carb versus SDC contrasts. | Both carbohydrate conditions showed higher SCFA trajectories than the no-added-carbohydrate condition in the descriptive data; adjusted pairwise contrasts will determine which condition-specific differences are statistically supported. |
| SDC produces the highest butyrate response. | Unsupported | The current mixed models test the three-level carbohydrate factor. The delta analysis reports an omnibus effect. Neither output establishes an adjusted SDC-versus-RDC butyrate difference. | Do not state until an adjusted SDC-versus-RDC contrast is reported with an effect estimate and confidence interval. |
| Obesity status does not substantially impair fermentation capacity. | Supported with qualification | No group main effect, group-by-carbohydrate interaction, group-by-time interaction, or three-way interaction remained significant after BH correction; the smallest adjusted group-related p-value was approximately 0.50 for butyrate. A nonsignificant test is not evidence of equivalence. | SCFA trajectories did not differ significantly by obesity status in this cohort. Avoid “equivalent,” “preserved,” or “unaffected” unless equivalence margins and confidence intervals are added. |
| Healthy participants had higher alpha diversity than participants with obesity at 48 hours under no added carbohydrate. | Unsupported after multiplicity correction | Raw Wilcoxon p-values were 0.0080 for Shannon, 0.0113 for Observed, and 0.0433 for Simpson, but BH-adjusted p-values were 0.102, 0.102, and 0.260, respectively. | No Healthy–Obese alpha-diversity comparison remained significant after FDR correction. |
| Carbohydrate, time, group, and group-by-carbohydrate each significantly affect beta diversity. | Unsupported as currently reported | The rendered PERMANOVA output provides only aggregate model results: all timepoints R²=0.171, p=0.001; 48 hours R²=0.0467, p=0.001. It does not expose term-specific R² or p-values. | Overall community composition differed across the modeled design. Term-specific effects require a complete PERMANOVA table. |
| Bifidobacterium increases with both carbohydrates, more under SDC. | Supported with qualification | MaAsLin2 coefficients versus the reference condition were 0.745 for SDC (q=0.023) and 0.582 for RDC (q=0.091). These are not a direct SDC-versus-RDC contrast. | *Bifidobacterium* was positively associated with carbohydrate exposure, with stronger evidence for SDC versus the reference condition. Do not claim SDC exceeds RDC without a direct contrast. |
| Differential abundance identifies multiple genera distinguishing RDC from SDC at 48 hours. | Unsupported in the plural | ANCOM-BC2 and MaAsLin2 each yielded one reported 48-hour RDC-versus-SDC genus at q<0.1. | One genus met the prespecified exploratory q<0.1 threshold in each 48-hour RDC-versus-SDC analysis. Report its identity, coefficient, interval, and method concordance. |
| Interpersonal variation persists after accounting for group and carbohydrate. | Supported with qualification | Subject-level trajectories and random-intercept models show substantial between-subject heterogeneity, but no variance decomposition or repeatability statistic is reported. | SCFA responses varied markedly among analyzed donor IDs. Quantify this with variance components or intraclass correlations before using “dominant individual axis.” |
| Butyrate and propionate responder pathways are partly independent. | Exploratory | The 2×2 overlap table is perfectly symmetric: 12 both responders, 12 both nonresponders, and 12 in each discordant cell. The denominator of 48 is subject-by-carbohydrate rows, not 48 independent participants, and the median split forces approximately equal groups within each stratum. | Median-defined butyrate and propionate response classifications showed incomplete overlap across subject-condition records. Do not infer pathway independence. |
| Butyrate responders have higher alpha diversity and distinct community structure. | Exploratory | Pooled RDC/SDC comparisons at 48 hours gave adjusted p=0.0167 for Observed, 0.00196 for Shannon, and 0.00348 for Simpson. Aggregate responder PERMANOVA gave R²=0.083, p=0.017. Responder labels and microbiome features come from the same experimental units, and carbohydrate-specific labels are pooled. | Median-defined butyrate response strata differed in 48-hour diversity in exploratory analyses. These analyses do not establish prediction and require carbohydrate-stratified, donor-aware validation. |
| Propionate responders have a weaker community signature. | Exploratory | Alpha metrics were not significant after adjustment (all adjusted p≥0.225); aggregate PERMANOVA gave R²=0.049, p=0.041. | Median-defined propionate strata showed no alpha-diversity difference and a small aggregate beta-diversity association. |
| Relative *Fusicatenibacter* abundance correlates with butyrate under SDC. | Unsupported as a significant association | At 48 hours under SDC, the linear model gave R²=0.0368, p=0.1569; Spearman ρ=0.287; N=56. Propionate results were R²≈0.004, p=0.655, ρ=-0.031. | Relative *Fusicatenibacter* abundance showed a weak, nonsignificant positive monotonic relationship with butyrate under SDC. |
| Absolute *Fusicatenibacter* abundance is higher in butyrate responders. | Exploratory and condition-dependent | The pooled RDC/SDC Wilcoxon comparison gave p=0.00123 with n=63 records per response group. SDC medians were 9.33×10⁶ copies/µL in responders and 0.978×10⁶ in nonresponders; RDC medians were 0.505×10⁶ and 0.196×10⁶. The p-value is pooled rather than SDC-specific, and multiple records per donor remain. | Absolute *Fusicatenibacter* abundance differed between median-defined butyrate response strata in a pooled exploratory analysis. A donor-aware response-by-carbohydrate model is required. |
| The *Fusicatenibacter* signal is pathway-specific. | Exploratory | The pooled propionate-responder comparison was nonsignificant (p=0.523; n=64 and 62 records). A significant-versus-nonsignificant comparison does not itself demonstrate different effects. | The absolute-abundance association was observed for butyrate-defined, but not propionate-defined, strata. Test a response-type interaction before claiming specificity. |
| Absolute *Bifidobacterium* abundance separates responders. | Unsupported | Pooled absolute-abundance tests were nonsignificant for butyrate (p=0.103) and propionate (p=0.531). Under SDC, the median was lower in butyrate responders (4.24×10⁷ copies/µL) than nonresponders (1.62×10⁸). | Absolute *Bifidobacterium* abundance did not differ significantly by median-defined responder status. |
| *Fusicatenibacter* is an SDC-utilizing butyrate producer. | Unsupported mechanistically | The study has no isolate growth, substrate-consumption, flux, metatranscriptomic, or co-culture data. The concurrent relative-abundance association with butyrate is nonsignificant. | The data nominate *Fusicatenibacter* for mechanistic follow-up; they do not establish SDC utilization or butyrate production by this genus. |
| The ex vivo assay and baseline microbiome can predict dietary response. | Unsupported | Response status is defined from the same ex vivo SCFA data used in downstream comparisons. Microbiome measurements are concurrent at 48 hours, not baseline predictors evaluated in an independent set. | The assay may support future hypothesis testing for response stratification. Prediction requires baseline-only features, held-out validation, and in vivo outcomes. |

## Supported numerical findings

### SCFA models

The integrated mixed model was:

`concentration ~ group * carbohydrate_type * timepoint_hr + (1 | subject)`

The strongest defensible statements are:

1. Time was associated with all five measured analytes after BH correction.
2. Carbohydrate condition was associated with acetate, butyrate, and propionate.
3. Carbohydrate-by-time interactions were supported for butyrate and propionate; the acetate interaction narrowly missed the adjusted threshold (adjusted p≈0.051).
4. No obesity-group main effect or group interaction was significant after BH correction.
5. The current model reports statistical significance but not the condition-specific estimated marginal means, effect sizes, or confidence intervals needed for a primary Results table.

### Community analyses

- Overall Bray–Curtis PERMANOVA, all timepoints: R²=0.171, p=0.001.
- Overall Bray–Curtis PERMANOVA, 48 hours: R²=0.0467, p=0.001.
- These are whole-model statistics. They cannot be assigned to obesity group, carbohydrate, time, or the group-by-carbohydrate interaction.
- Healthy–Obese alpha-diversity comparisons did not survive BH correction.
- RDC-versus-SDC alpha-diversity comparisons at 48 hours were nonsignificant within both groups.

### Responder analyses

- “Responder” means ΔSCFA at or above the median within an analyte-by-carbohydrate stratum.
- The classification is relative and sample-dependent. It is not a biological or clinical threshold.
- The overlap denominator reported as “48 subjects” is 48 subject-by-carbohydrate records.
- ANCOM-BC2 identified zero responder-associated taxa at q<0.1 for both butyrate and propionate; MaAsLin2 reported one association for each.
- All responder analyses should be labeled secondary and hypothesis-generating.

## Methods and sample-accounting audit

### Planned and analyzed units

| Level | Current evidence | Required interpretation |
|---|---|---|
| Planned cohort | 40 participants: 20 healthy-weight and 20 with obesity | Protocol target, not analyzed N |
| SCFA source labels | 32 A/B-coded labels before collapse to 16 numeric subject IDs: eight control and eight case | The biological meaning of A/B is undocumented |
| Microbiome observations | 292 sample-level profiles in the phyloseq object, all matched to qPCR data | Independent donor, collection, and well-replicate counts must be resolved |
| Responder overlap | 48 subject-by-carbohydrate records | Must not be described as 48 participants |

The SCFA and integrated workflows define subject as:

`str_extract(sampleid, "^[0-9]+")`

This removes the A/B suffix before averaging technical replicates and fitting random-intercept models. The manuscript cannot report an analyzed participant count until the study team confirms whether A/B denotes independent donors, repeated collections, split aliquots, or another unit.

### Timepoint conflict

The integrated Methods state that SCFA was measured at 0, 48, and 72 hours. The analyzed files and primary SCFA workflow contain 0- and 48-hour records. No 72-hour result is present in the integrated analysis. The manuscript should use 0 and 48 hours unless a separate 72-hour dataset is supplied and analyzed.

### Terminology conflicts

- The project description uses “fast digestible carbohydrate” (FDC), whereas the executed analysis uses “rapid digestible carbohydrate” (RDC). Confirm that these refer to the same substrate before standardizing on RDC.
- The project description mentions untreated and prebiotic GFA-treated samples; the executed analysis uses no added carbohydrate, RDC, and SDC. Confirm whether GFA was omitted, renamed, or represented by one of the analyzed conditions.
- Analysis labels alternate between Control/Case and Healthy/Obese. The manuscript should use healthy-weight and obesity groups after the cohort definitions and mapping are verified.
- Use 0 hours and 48 hours consistently; do not retain the unobserved 72-hour label.

### Missing cohort information

The repository does not establish:

- analyzed age range, mean age, or age distribution;
- BMI percentile or other obesity case definition;
- inclusion and exclusion criteria;
- sex, race/ethnicity, medications, recent antibiotic exposure, or relevant diet variables;
- recruitment dates and setting details;
- ethics board, protocol number, consent, and assent statements.

### Missing fermentation protocol

The repository does not establish:

- exact RDC and SDC chemical identities, source, lot, purity, or concentration;
- culture-medium composition;
- inoculum preparation and fecal dilution;
- vessel or well volume;
- anaerobic atmosphere, incubation temperature, shaking, and pH control;
- treatment allocation, plate design, batch handling, and the rationale for r1/r2 and s1/s2;
- whether “no added carbohydrate” retained background carbohydrate in the medium.

### SCFA assay

Locally documented elements include PFBBr derivatization, GC-negative chemical ionization-MS, quantitative acetate/propionate/butyrate calibration, and a succinate calibration standard. The detailed facility method specifies 100 µL extract, 100 mM borate buffer, 100 mM PFBBr, hexane extraction, 65°C for one hour, and MassHunter B.10.

Items requiring reconciliation:

- The overview names an Agilent 8890; the detailed method names a 7890A GC with a 5975C detector.
- The facility SOP describes extraction from 100 mg fecal/cecal material per mL of 80% methanol, while the submission sheet identifies 500 mg fecal-culture material.
- The manuscript must distinguish quantitatively calibrated SCFAs from compounds reported by normalized abundance.
- The analysis reports concentrations in µM; this unit should be verified against the exported quantification file and facility report.

### 16S and absolute abundance

The integrated workflow uses Zymo Bac16S V3–V4 data, DADA2-derived ASVs, genus-level aggregation, phyloseq, ANCOM-BC2, and MaAsLin2. The exact sequencing primers, DADA2 filtering parameters, taxonomy database/version, library-depth exclusions, contamination handling, and run accession are not yet manuscript-ready.

Absolute taxon abundance is calculated as relative abundance multiplied by qPCR-derived total bacterial genome copies per µL. The manuscript should report copies per µL, not “copies per sample,” and add the qPCR target, primer sequences, efficiency, standard curve, limit of quantification, extraction input, and the gene-to-genome conversion assumptions.

## Analysis gate before claim-locked drafting

### Priority 0: blocking

1. **Resolve biological units.** Document the meaning of A/B and rebuild donor, collection, aliquot, and technical-replicate identifiers.
2. **Recompute analyzed N.** Provide a participant flow and condition-by-timepoint completeness table after exclusions.
3. **Obtain the fermentation protocol.** Exact substrate composition and culture conditions are necessary to interpret “digestibility.”
4. **Run direct primary contrasts.** Report adjusted RDC versus SDC, RDC versus no-carb, and SDC versus no-carb estimates with 95% confidence intervals for each primary SCFA.
5. **Export term-level PERMANOVA.** Report R² and p for each model term and test multivariate dispersion.
6. **Correct the Fusicatenibacter summary.** Remove “significant correlation” language for the relative-abundance analysis.

### Priority 1: required for responder and focal-taxon claims

1. Model SCFA deltas with donor-aware repeated measures rather than a combined one-way ANOVA.
2. Treat SCFA response continuously as the primary interpersonal-variation analysis.
3. If median groups are retained, stratify or model interaction by carbohydrate and account for repeated records per donor.
4. Compare response-type effects directly before claiming butyrate specificity.
5. Use baseline microbiome features for any prediction analysis and validate out of sample.
6. Prespecify the multiplicity families for SCFA, diversity, focal taxa, and discovery-wide taxon tests.

### Priority 2: strengthening

1. Report variance components or intraclass correlations for donor heterogeneity.
2. Add sensitivity analyses using one microbiome profile per donor-condition-time unit.
3. State which findings reproduce across ANCOM-BC2 and MaAsLin2.
4. Report effect estimates and confidence intervals alongside p/q values.
5. If feasible, add strain-level substrate utilization or co-culture validation before making a mechanistic claim.

## Language guardrails

Use:

- “associated with”
- “did not differ significantly”
- “median-defined response strata”
- “concurrent 48-hour association”
- “hypothesis-generating”
- “nominates for follow-up”

Avoid until supported:

- “governs”
- “supersedes obesity status”
- “preserved capacity” as an equivalence claim
- “defines responder phenotypes”
- “predicts response”
- “pathway-specific”
- “SDC-utilizing butyrate producer”
- “mechanistic keystone”
- “companion diagnostic”

## Recommended title after audit

**Preferred current title**

> Carbohydrate Condition and Incubation Time Shape Ex Vivo Short-Chain Fatty Acid Production by Adolescent Fecal Microbiota with Limited Modification by Obesity Status

**Alternative if direct RDC-versus-SDC contrasts support digestibility-specific effects**

> Digestible Carbohydrate Type Shapes Ex Vivo Short-Chain Fatty Acid Production and Microbiome Responses in Adolescents with and without Obesity

Do not place *Fusicatenibacter*, “responder phenotypes,” or precision nutrition in the title before donor-aware, condition-specific, and independent validation.

## Source ledger

- `integrated/integrated-scfa-microbiome.Rmd`: sample-ID parsing and averaging, lines 126–184; Methods, lines 350–415; sample accounting, lines 418–485; SCFA mixed models, lines 749–812; beta-diversity models, lines 1133–1362; Fusicatenibacter correlations, lines 4322–4453; responder definition and overlap, lines 4480–4570; responder community models, lines 4573–4765 and 5180–5330; absolute-abundance responder analyses, lines 5716 onward.
- `integrated/integrated-scfa-microbiome.html`: rendered mixed-model table in section 5.4; diversity output in sections 6.2–6.3; Fusicatenibacter correlations in section 7.4; responder output in sections 7.5–7.9.
- `scfa_metabolomics/scfa-project2-analysis-improved-v2.Rmd`: sample-ID parsing and technical-replicate averaging, lines 174–229.
- `docs/2025_DFI_HMMF_Methods.pdf`: PFBBr panel and detailed GC-negative chemical ionization-MS method, pages 2–3.
- `docs/Introduction-information.md`: planned cohort and aims, lines 158–181.
- `docs/zc07-next-steps-proposal.2026-02-01.md`: current narrative claims, limitations, and manuscript proposal.

## Go/no-go checklist

### Begin full manuscript drafting when

- [ ] A/B and technical-replicate structure is documented.
- [ ] Analyzed participant N and exclusions are reconciled.
- [ ] Direct adjusted carbohydrate contrasts are available.
- [ ] Term-level PERMANOVA and dispersion results are available.
- [ ] Exact fermentation and substrate methods are supplied.
- [ ] Cohort and ethics fields are supplied.
- [ ] Responder language is explicitly exploratory or independently validated.
- [ ] Focal-taxon claims use condition-specific, donor-aware models.

### Position for Gut Microbes when

- [ ] The primary story rests on supported SCFA and community effects rather than a single exploratory taxon.
- [ ] Absolute-abundance methods are fully reproducible.
- [ ] Mechanistic language is removed or supported by strain/co-culture data.
- [ ] All principal results include effect sizes, confidence intervals, and corrected p/q values.
- [ ] The Discussion clearly separates ex vivo metabolic capacity from in vivo response.

**Current status:** proceed with the manuscript framework and targeted literature work; do not lock the title, Abstract, or focal-taxon conclusion until Priority 0 items are resolved.
