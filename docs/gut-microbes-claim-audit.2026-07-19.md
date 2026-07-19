# Gut Microbes manuscript claim audit

**Project:** Abbott–Lurie ZC07, Project 2  
**Date:** July 19, 2026  
**Purpose:** Determine which proposed manuscript claims are supported by the current analysis, which require qualification, and which require reanalysis or additional source information before drafting a *Gut Microbes* submission.

## Executive decision

The current dataset supports an ex vivo phenotyping paper centered on incubation time, carbohydrate condition, limited modification of SCFA trajectories by obesity status, and marked interpersonal variation. It does not yet support the proposed mechanistic title or a validated responder-prediction claim.

The manuscript can proceed to a structured shell, literature mapping, and provisional Methods. Claim-locked Results and Abstract drafting should wait until facility unit/QC confirmation and the primary RDC-versus-SDC contrasts are refreshed after the nested-unit model.

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
| SDC produces the highest butyrate response. | Unsupported | With nested random effects for subject, A/B aliquot, and culture well (re-knit 2026-07-19), the direct mixed-model contrast estimated the SDC-versus-RDC difference in 0-to-48-hour butyrate change as 1.03 µM (95% CI, −0.85 to 2.91; Tukey-adjusted p=0.396). The concentration unit remains pending facility confirmation. | The butyrate increase was numerically greater under SDC than RDC, but the adjusted contrast was not significant and does not demonstrate SDC superiority. |
| Obesity status does not substantially impair fermentation capacity. | Supported with qualification | After clustering on 16 independent numeric subjects (8 per group), with A/B aliquots and culture wells nested within subject, none of the six primary obesity-minus-healthy-weight differences in 0-to-48-hour change was significant. Estimates ranged from −1.19 to 3.43 in the concentration units labeled µM. No external biological equivalence margin or project-specific assay-precision bound was available, so formal TOST equivalence was not evaluable. | SCFA trajectories did not differ significantly by obesity status in this cohort (N=16 subjects). Equivalence was not evaluated because no externally justified margin was available; avoid “equivalent,” “preserved,” or “unaffected.” |
| Healthy participants had higher alpha diversity than participants with obesity at 48 hours under no added carbohydrate. | Unsupported after multiplicity correction | After donor-level aggregation (one mean per donor × carbohydrate × timepoint), raw Wilcoxon p-values at 48 h under no added carbohydrate were 0.0104 for Shannon, 0.0070 for Observed, and 0.105 for Simpson; BH-adjusted values were approximately 0.094, 0.094, and 0.47. | No Healthy–Obese alpha-diversity comparison remained significant after FDR correction. |
| Carbohydrate, time, group, and group-by-carbohydrate each significantly affect beta diversity. | Unsupported as currently reported | Donor-aware PERMANOVA (`strata = donor_id`) gave whole-model summaries: all timepoints R²=0.174, p=0.001; 48 hours R²=0.053, p=0.004. Term-specific R²/p are still not exported in the rendered table. | Overall community composition differed across the modeled design under donor-level strata. Term-specific effects still require a complete PERMANOVA term table. |
| Bifidobacterium increases with both carbohydrates, more under SDC. | Supported with qualification | MaAsLin2 coefficients versus the reference condition were 0.745 for SDC (q=0.023) and 0.582 for RDC (q=0.091). These are not a direct SDC-versus-RDC contrast. | *Bifidobacterium* was positively associated with carbohydrate exposure, with stronger evidence for SDC versus the reference condition. Do not claim SDC exceeds RDC without a direct contrast. |
| Differential abundance identifies multiple genera distinguishing RDC from SDC at 48 hours. | Unsupported in the plural | ANCOM-BC2 and MaAsLin2 each yielded one reported 48-hour RDC-versus-SDC genus at q<0.1. | One genus met the prespecified exploratory q<0.1 threshold in each 48-hour RDC-versus-SDC analysis. Report its identity, coefficient, interval, and method concordance. |
| Interpersonal variation persists after accounting for group and carbohydrate. | Supported with qualification | Subject-level trajectories and random-intercept models show substantial between-subject heterogeneity, but no variance decomposition or repeatability statistic is reported. | SCFA responses varied markedly among analyzed donor IDs. Quantify this with variance components or intraclass correlations before using “dominant individual axis.” |
| Butyrate and propionate responder pathways are partly independent. | Exploratory | The 2×2 overlap table remains symmetric after the donor-level rebuild: 12 both responders, 12 both nonresponders, and 12 in each discordant cell. The denominator is **48 donor × carbohydrate records** (16 donors × 3 carbohydrates), not 48 people, and the median split forces approximately equal groups within each stratum. | Median-defined butyrate and propionate response classifications showed incomplete overlap across donor-condition records. Do not infer pathway independence. |
| Butyrate responders have higher alpha diversity and distinct community structure. | Exploratory | Pooled RDC/SDC comparisons at 48 hours gave adjusted p=0.0167 for Observed, 0.00196 for Shannon, and 0.00348 for Simpson. Aggregate responder PERMANOVA gave R²=0.083, p=0.017. Responder labels and microbiome features come from the same experimental units, and carbohydrate-specific labels are pooled. | Median-defined butyrate response strata differed in 48-hour diversity in exploratory analyses. These analyses do not establish prediction and require carbohydrate-stratified, donor-aware validation. |
| Propionate responders have a weaker community signature. | Exploratory | Alpha metrics were not significant after adjustment (all adjusted p≥0.225); aggregate PERMANOVA gave R²=0.049, p=0.041. | Median-defined propionate strata showed no alpha-diversity difference and a small aggregate beta-diversity association. |
| Relative *Fusicatenibacter* abundance correlates with butyrate under SDC. | Unsupported as a significant association | At 48 hours under SDC, the linear model gave R²=0.0368, p=0.1569; Spearman ρ=0.287; N=56. Propionate results were R²≈0.004, p=0.655, ρ=-0.031. | Relative *Fusicatenibacter* abundance showed a weak, nonsignificant positive monotonic relationship with butyrate under SDC. |
| Absolute *Fusicatenibacter* abundance is higher in butyrate responders. | Exploratory and condition-dependent | The pooled RDC/SDC Wilcoxon comparison gave p=0.00123 with n=63 records per response group. SDC medians were 9.33×10⁶ copies/µL in responders and 0.978×10⁶ in nonresponders; RDC medians were 0.505×10⁶ and 0.196×10⁶. The p-value is pooled rather than SDC-specific, and multiple records per donor remain. | Absolute *Fusicatenibacter* abundance differed between median-defined butyrate response strata in a pooled exploratory analysis. A donor-aware response-by-carbohydrate model is required. |
| qPCR-scaled *Fusicatenibacter* net change associates with SDC butyrate net change. | Unsupported as significant; exploratory | Donor-level SDC-minus-no-carb Δlog10(qPCR-scaled *Fusicatenibacter*) vs Δbutyrate: Spearman ρ≈0.48, permutation p≈0.088, 95% bootstrap CI −0.24 to 0.89, N=14 donors (`integrated/results/additional_analyses_2026-07-19/focal_fusicatenibacter/`). Specificity family BH q-values all ≥0.95. | Directionally positive concurrent association after correct nesting; interval includes zero. Label exploratory and say “qPCR-scaled,” not genus-specific qPCR. |
| The *Fusicatenibacter* signal is pathway-specific. | Exploratory | The pooled propionate-responder comparison was nonsignificant (p=0.523; n=64 and 62 records). A significant-versus-nonsignificant comparison does not itself demonstrate different effects. | The absolute-abundance association was observed for butyrate-defined, but not propionate-defined, strata. Test a response-type interaction before claiming specificity. |
| Absolute *Bifidobacterium* abundance separates responders. | Unsupported | Pooled absolute-abundance tests were nonsignificant for butyrate (p=0.103) and propionate (p=0.531). Under SDC, the median was lower in butyrate responders (4.24×10⁷ copies/µL) than nonresponders (1.62×10⁸). | Absolute *Bifidobacterium* abundance did not differ significantly by median-defined responder status. |
| *Fusicatenibacter* is an SDC-utilizing butyrate producer. | Unsupported mechanistically | The study has no isolate growth, substrate-consumption, flux, metatranscriptomic, or co-culture data. The concurrent relative-abundance association with butyrate is nonsignificant. | The data nominate *Fusicatenibacter* for mechanistic follow-up; they do not establish SDC utilization or butyrate production by this genus. |
| The ex vivo assay and baseline microbiome can predict dietary response. | Unsupported | Response status is defined from the same ex vivo SCFA data used in downstream comparisons. Microbiome measurements are concurrent at 48 hours, not baseline predictors evaluated in an independent set. | The assay may support future hypothesis testing for response stratification. Prediction requires baseline-only features, held-out validation, and in vivo outcomes. |

## Supported numerical findings

### SCFA models

The corrected integrated mixed model is:

`concentration ~ group * carbohydrate_type * timepoint_hr + (1 | subject) + (1 | aliquot_id) + (1 | well_id)`

where `subject` is the numeric donor ID, `aliquot_id` is the A/B sample split, and `well_id` is the culture well within aliquot × carbohydrate.

The strongest defensible statements are:

1. The previous models that treated A/B labels as independent subjects overstated the independent N; clustering must be on the numeric subject.
2. The corrected model analyzes 16 independent subjects, 32 A/B aliquots, and 160 culture wells.
3. The six primary obesity-group difference-in-change estimates and confidence intervals are available.
4. None of those six superiority contrasts was significant; this does not establish equivalence.

### Obesity-group equivalence and non-inferiority analysis

The corrected analysis defines the estimand separately within each carbohydrate condition as:

`(Obesity 48 h − Obesity 0 h) − (Healthy-weight 48 h − Healthy-weight 0 h)`

The nested model uses subject → A/B aliquot → culture well. The six primary estimates, with model-based 95% confidence intervals, were:

- acetate RDC: −0.09 (−5.17 to 5.00);
- acetate SDC: 3.43 (−1.66 to 8.51);
- propionate RDC: 0.74 (−0.70 to 2.18);
- propionate SDC: 0.38 (−1.06 to 1.82);
- butyrate RDC: −0.76 (−2.96 to 1.44);
- butyrate SDC: −1.19 (−3.39 to 1.02).

Values use the concentration units labeled µM in the current workflow; facility confirmation is still required. Donor-level bootstrap intervals (resampling N=16 subjects) gave the same qualitative conclusion. Model diagnostics showed non-normal residuals for all three primary analytes and heteroscedasticity for propionate and butyrate, so both model-based and bootstrap intervals should be retained.

According to PubMed, available PFBBr assay papers provide analytical detection, recovery, or precision information but not a biologically negligible difference for this ex vivo obesity-group estimand ([DOI](https://doi.org/10.1016/j.jchromb.2018.06.028); [DOI](https://doi.org/10.1016/j.jchromb.2023.123826)). Comparable obesity fermentation studies also did not supply a reusable equivalence margin ([DOI](https://doi.org/10.1128/mBio.00914-20); [DOI](https://doi.org/10.3390/nu11020217)).

The primary margin file is therefore deliberately blank. TOST p-values and pass/fail decisions are not reported. The post hoc symmetric bounds required merely to contain the observed 90% intervals ranged from ±1.59 to ±7.70 µM across the six contrasts; these are precision frontiers, not margins to adopt after seeing the data.

**Decision:** Retain “did not differ significantly” and state that equivalence was not evaluable because an externally justified margin was unavailable. Do not use “preserved,” “equivalent,” “non-inferior,” or “unaffected.”

### Community analyses

- Overall Bray–Curtis PERMANOVA (donor-level strata), all timepoints: R²=0.174, p=0.001.
- Overall Bray–Curtis PERMANOVA (donor-level strata), 48 hours: R²=0.053, p=0.004.
- These remain whole-model statistics. They cannot be assigned to obesity group, carbohydrate, time, or the group-by-carbohydrate interaction until term-level tables are exported.
- Healthy–Obese alpha-diversity comparisons after donor aggregation did not survive BH correction.
- RDC-versus-SDC alpha-diversity comparisons at 48 hours should be re-checked from the donor-aware re-knit before quoting.

### Responder analyses

- “Responder” means ΔSCFA at or above the median within an analyte-by-carbohydrate stratum at the **donor** grain.
- The classification is relative and sample-dependent. It is not a biological or clinical threshold.
- The overlap denominator is 48 donor × carbohydrate records (16 × 3).
- ANCOM-BC2 / MaAsLin2 responder associations from the donor-aware re-knit remain exploratory; confirm method concordance before citing taxa.
- All responder analyses should be labeled secondary and hypothesis-generating.

## Methods and sample-accounting audit

### Planned and analyzed units

| Level | Current evidence | Required interpretation |
|---|---|---|
| Planned cohort | 20 participants in submission metadata (40 A/B aliquot labels when each person has A and B) | Protocol target, not analyzed N |
| SCFA subjects | 16 analyzed numeric donors: 8 healthy-weight and 8 obesity | Independent biological unit for group inference |
| SCFA aliquots | 32 A/B sample splits (2 per analyzed subject) | Nested within subject; not independent people |
| SCFA culture wells | R1/R2 and S1/S2 wells plus one no-added-carbohydrate well per aliquot-condition | Nested within aliquot × carbohydrate |
| Microbiome observations | 292 sample-level profiles in the phyloseq object, all matched to qPCR data | Match on A/B aliquot labels; cluster inference by numeric subject |
| Responder overlap | aliquot-by-carbohydrate records | Must not be described as independent participants |

The corrected SCFA workflow defines:

`subject = str_extract(sampleid, "^[0-9]+")`  
`aliquot = str_to_upper(str_match(sampleid, "^[0-9]+([A-Za-z])")[, 2])`

The facility metadata contain 40 A/B aliquot labels from 20 numeric subjects. The SCFA export contains 32 aliquots from 16 subjects. The remaining accounting task is to document why four planned subjects (eight A/B labels) are absent from the SCFA export and why two obesity no-added-carbohydrate 48-hour observations are missing.

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

> **Unit-hierarchy rebuild (2026-07-19):** Baseline preserved at git `67bd320`. Donor-aware re-knit completed (`render/integrated-scfa-microbiome.html`). Community/responder cells below were refreshed from that render where numbers are quoted; remaining exploratory DA taxa lists and Fusicatenibacter responder p-values should be treated as provisional until method-by-method concordance is reconfirmed in a dedicated pass. Nested Python SCFA obesity contrasts (N=16) remain the SCFA group authority. See [`docs/unit-hierarchy-reanalysis-plan.2026-07-19.md`](unit-hierarchy-reanalysis-plan.2026-07-19.md) and `integrated/metadata/`.
>
> **Additional analyses execution (2026-07-19):** Core Priority A/B analyses from [`docs/additional-analyses-plan.2026-07-19.md`](additional-analyses-plan.2026-07-19.md) were run into `integrated/results/additional_analyses_2026-07-19/`. Highlights: adjusted carbohydrate-change contrasts for acetate/propionate/butyrate; donor-aware alpha models; term-level PERMANOVA + betadisper on donor-aggregated profiles; exact well-matched qPCR-scaled *Fusicatenibacter* donor-level net-change analysis (ρ≈0.48, permutation p≈0.088, bootstrap CI includes 0, N=14). MaAsLin 3 / ALDEx3 await a separate R≥4.6 environment. Facility qPCR LOD/LOQ docs remain outstanding. Do not claim a significant *Fusicatenibacter*–butyrate association from the new donor-level test.

### Priority 0: blocking

1. **Reconcile exclusions.** The corrected SCFA analysis contains 16 independent subjects (8 per group) and 32 A/B aliquots from a planned 20-subject submission sheet; document why four planned subjects are absent and why two obesity no-added-carbohydrate 48-hour observations are missing.
2. **Confirm assay units and QC.** Obtain facility confirmation of the exported concentration unit, dilution mapping, analyte-specific CV, and LOQ.
3. **Obtain the fermentation protocol.** Exact substrate composition and culture conditions are necessary to interpret “digestibility.”
4. **Run direct primary contrasts.** Report adjusted RDC versus SDC, RDC versus no-carb, and SDC versus no-carb estimates with 95% confidence intervals for each primary SCFA.
5. **Export term-level PERMANOVA.** Report R² and p for each model term and test multivariate dispersion.
6. **Correct the Fusicatenibacter summary.** Remove “significant correlation” language for the relative-abundance analysis.
7. **Lock equivalence margins only if externally justified.** Until study-team biological thresholds and project QC are available, treat equivalence as not evaluable.

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
- “equivalence was not evaluable because no external margin was available”
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

- `integrated/integrated-scfa-microbiome.Rmd`: nested SCFA + microbiome `donor_id` clustering, donor-aggregated alpha tests, donor-aware RE/strata.
- `render/integrated-scfa-microbiome.html` (and copy under `integrated/`): 2026-07-19 donor-aware re-knit; prefer this over any earlier HTML/PDF.
- `integrated/metadata/`: Phase 0 inventory and canonical experimental-unit tables.
- `scfa_metabolomics/scfa-project2-analysis-improved-v2.Rmd`: legacy standalone workflow that collapses A/B and averages R/S wells; retained as historical source only.
- `scfa_metabolomics/obesity_equivalence_analysis.py`: corrected nested-well trajectory models, 90% and 95% model intervals, donor-level bootstrap sensitivity intervals, diagnostics, influence checks, and conditional TOST calculations.
- `scfa_metabolomics/equivalence_margins.csv`: machine-readable margin lock; version `v0.1-unlocked` contains no numerical bounds.
- `scfa_metabolomics/results/obesity_group_scfa_*.csv`: complete contrasts, diagnostics, sample counts, and analysis status.
- `docs/scfa-obesity-equivalence-margin-spec.2026-07-19.md`: unit review, external-evidence review, margin decision, and requirements for a future locked version.
- `docs/2025_DFI_HMMF_Methods.pdf`: PFBBr panel and detailed GC-negative chemical ionization-MS method, pages 2–3.
- `docs/Introduction-information.md`: planned cohort and aims, lines 158–181.
- `docs/zc07-next-steps-proposal.2026-02-01.md`: current narrative claims, limitations, and manuscript proposal.

## Go/no-go checklist

### Begin full manuscript drafting when

- [x] Subject → A/B aliquot → culture-well nesting is documented and implemented for SCFA.
- [x] Microbiome `donor_id` demotes A/B from the top-level clustering unit; RE/strata use `donor_id`.
- [x] Donor-aware integrated HTML re-knit produced (`render/integrated-scfa-microbiome.html`).
- [ ] Equivalence margins are locked from external evidence (still unavailable).
- [ ] Facility unit / CV / LOQ confirmation is complete.
- [ ] Term-level PERMANOVA table is exported (still whole-model only).
- [ ] Responder language is explicitly exploratory or independently validated.
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
