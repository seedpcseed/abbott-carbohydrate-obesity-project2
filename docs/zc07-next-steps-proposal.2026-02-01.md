---
title: "Abbott-Lurie Children's Collaboration ZC07"
subtitle: "Next Steps Proposal"
date: "February 1, 2026"
mainfont: "Roboto"
geometry: "top=1in, bottom=1in, left=0.5in, right=0.5in"
output: pdf_document
header-includes:
  - \usepackage{fontspec}
  - \setmainfont{DejaVu Sans}
  - \defaultfontfeatures{Ligatures=TeX}
---

# Abbott-Lurie Children's Collaboration ZC07 

## Introduction and Purpose 

This document provides a comprehensive review of the Abbott-Lurie Children's Collaboration ZC07 project, summarizing progress to date and outlining a strategic roadmap for near-term, mid-term, and long-term activities. The objectives of this proposal are to:

1. **Review project goals and milestones:** Summarize the original aims and assess progress toward each objective.
2. **Identify completed and ongoing work:** Document analyses and experiments that have been finalized, and clarify remaining tasks.
3. **Address challenges:** Describe barriers encountered in completing outstanding aims and propose solutions.
4. **Outline dissemination strategy:** Present a structured abstract for conference presentation and a publication outline with key content, figures, and tables.
5. **Define mechanistic follow-on studies:** Identify experiments required to strengthen publication impact by elucidating mechanisms underlying observed phenotypes.
6. **Propose mid- and long-term investigations:** Articulate a research trajectory that extends current findings toward translational and clinical applications. 

## Project Background and Progress to Date 

**Project goals:**

- Measure variation in responses of human gut-associated microbiome communities to fast (FDC/RDC) and slow digestible (SDC) carbohydrates.
- Identify childhood-associated organisms with facile utilization of SDC.
- Test interindividual variation of short-chain fatty acid (SCFA) production among fecal microbiota samples to slow and fast digestible carbohydrates.
- Measure energy harvest differences between human obesity-associated fecal microbiota (Aim 2).
- Isolate SDC bacterial utilizers using single-cell isolation toward synbiotic combinations (Aim 3).

**Microbiome analysis progress (microbiome_analysis/analysis.v4_zymo.html):**

- 16S rRNA gene sequencing (Zymo pipeline) with phyloseq-based analysis; alpha diversity (Observed, Shannon, Simpson) at 0H and 48H; beta diversity (Bray-Curtis PCoA, PERMANOVA with subject strata); baseline clustering and trajectories; genus-level composition changes 0H->48H by carbohydrate type; expanded differential abundance (corncob, ANCOM-BC2, LEfSe, MaAsLin2, etc.) with 48H focus; FDR control throughout.

**SCFA analysis progress (scfa_metabolomics/scfa-project2-analysis-improved-v2.html):**

- Targeted SCFA metabolomics (GC-MS/PFBBr: acetate, butyrate, propionate, 5-aminovalerate, succinate); summary statistics by group and carbohydrate type; SCFA concentrations over time; delta (48h-0h) response magnitude; two-way ANOVA (group × carbohydrate) on deltas; carbohydrate type effects (combined and stratified by Control/Case); mixed-effects models for SCFA trajectories. Key outcomes: carbohydrate type and time drive SCFA production; RDC and SDC both increase SCFAs vs no-carbohydrate; limited Case vs Control differences; pronounced interpersonal variation in SCFA response.

**Integrated analysis (integrated/integrated-scfa-microbiome.html):**

- Single workflow combining SCFA and microbiome with shared metadata; SCFA delta and mixed-effects results; alpha/beta diversity; composition and differential abundance; taxon–SCFA associations (e.g., *Fusicatenibacter* vs butyrate/propionate); responder vs non-responder (median-delta split for butyrate/propionate) and composition comparisons; absolute abundance (qPCR-based) for focal taxa; correlation and regression for taxa–SCFA relationships.

## Steps in Progress

- **Aim 2 (energy harvest):** Ex vivo exposure of fecal samples to SDC and FDC with timepoint collection and bomb calorimetry (kcal per gram) not yet completed.
- **Aim 3 (SDC utilizers):** Single-cell capture, arraying, and growth assays for SDC vs FDC utilizers from high-SCFA / low–energy-extraction donors not yet completed.
    - We are currently obtaining fully sequenced reference strains of *Fusicatenibacter* (details in `integrated/fusicatenibacter_strain_information_DSMZ.md`). Two type strains are available: *F. saccharivorans* (DSM 26062, HT03-11), isolated from human feces, with documented saccharolytic activity (glucose, arabinose, xylose, cellobiose, starch) and relevance to carbohydrate metabolism and SCFA production in Lachnospiraceae; and *F. faecihominis* (DSM 113288), also from human feces, with genome sequence available for comparative and functional validation.
    - We may want to obtain key species of *Bifidobacterium* as well, since abundance changed with both RDC and SDC (more so with SDC). *B. adolescentis* and *B. longum* are well-documented cross-feeders for butyrate producers: they ferment oligosaccharides and starch to lactate and acetate, which lactate-utilizing taxa (e.g., *Anaerostipes*, *Eubacterium hallii*, *Roseburia*, *Faecalibacterium prausnitzii*) convert to butyrate; they can also release partially degraded oligosaccharides used directly by some butyrate producers. Reference strains of these species would support mechanistic and synbiotic work alongside *Fusicatenibacter*.
- **Analysis:** Finalizing responder/non-responder definitions and validation; any remaining cross-dataset QC and narrative for publication.

**Challenges in completing remaining goals:** Resource and timeline constraints for bomb calorimetry and single-cell isolation workflows; need to align sample availability and lab capacity for Aims 2–3; defining responder thresholds (median-delta vs clinical anchors) for publication; coordinating integrated analyses across SCFA and microbiome pipelines.

## Near Term Work 
### Abstract (Sample)
**Introduction.** Childhood obesity is a major driver of adult obesity, and the gut microbiome is a central mediator of metabolic health. Carbohydrate quality—particularly digestibility—may matter more than quantity for metabolic outcomes. Interpersonal variation in microbial short-chain fatty acid (SCFA) production may determine who benefits from slowly digestible carbohydrates (SDC). We aimed to characterize ex vivo SCFA production and microbiome composition in adolescents with and without obesity in response to different carbohydrate substrates.

**Methods.** Fecal samples from adolescents with obesity (Case) and healthy-weight controls (Control) were incubated ex vivo under no-carbohydrate, rapid digestible (RDC), and slow digestible (SDC) conditions. SCFAs (acetate, butyrate, propionate, and related metabolites) were quantified by gas chromatography–mass spectrometry after PFBBr derivatization; microbiome composition was assessed by 16S rRNA gene sequencing. Analyses included linear mixed-effects models for SCFA trajectories, alpha and beta diversity with PERMANOVA, differential abundance (ANCOM-BC2, MaAsLin2), taxon–SCFA associations, and responder vs non-responder stratification (median-delta SCFA).

**Results.** Carbohydrate type and incubation time were the primary drivers of SCFA production; both RDC and SDC increased SCFAs versus no-carbohydrate controls. Obesity status did not markedly blunt fermentation capacity. Pronounced interpersonal variation in SCFA response was observed. Specific taxa (e.g., *Fusicatenibacter*) were associated with SCFA concentrations; responder stratification enabled targeted characterization of microbiome composition differences.

**Conclusions.** Ex vivo fermentation capacity for RDC and SDC is preserved in adolescents with obesity, with substantial interindividual variation in SCFA production. Taxon–SCFA associations and responder phenotypes support precision nutrition strategies and future trials stratified by microbiome/SCFA phenotype.

**Target audience/meetings:** Pediatric obesity, nutrition, microbiome (e.g., DDW, PAS, Obesity Week, AGA, ASM Microbe); Abbott-internal and Lurie–Abbott collaboration meetings.

### Publication 

- **Title/structured abstract:** Short, descriptive title; structured abstract (Background, Methods, Results, Conclusions).
- **Introduction:** Childhood obesity, microbiome as metabolic mediator, SDC/FDC rationale, Lurie–Abbott aims, interpersonal variation and precision nutrition.
- **Methods:** Subject/sample description; ex vivo culture conditions; SCFA (extraction, derivatization, GC-MS); 16S sequencing and processing; alpha/beta diversity, differential abundance, mixed-effects and taxon–SCFA analyses; responder definition; software and FDR.
- **Results (content we have):** Sample summary; SCFA by group and carbohydrate type; SCFA over time and delta analyses; mixed-effects results; alpha/beta diversity; composition and differential abundance; taxon–SCFA (e.g., *Fusicatenibacter*) and responder vs non-responder composition.
- **Results (to fill in):** Bomb calorimetry (if Aim 2 completed); isolate growth data (if Aim 3 used); any additional validation or sensitivity analyses.
- **Discussion:** Interpretation of carbohydrate and time effects; interpersonal variation and precision nutrition; taxa–SCFA and responder findings; limitations (ex vivo, single cohort, timepoints); future directions (Aim 2/3, mechanisms, trials).
- **Major figures/tables planned:** (1) Study design/sampling schema; (2) SCFA by carbohydrate type and time (e.g., line + delta boxplots); (3) Alpha/beta diversity by group, carbohydrate, timepoint; (4) Composition (e.g., genus bar/alluvial) 0H vs 48H; (5) Differential abundance (forest or volcano); (6) Taxon–SCFA correlations (e.g., *Fusicatenibacter* vs butyrate/propionate); (7) Responder vs non-responder composition; (8) Summary tables: mixed-effects SCFA, PERMANOVA, DA summary (taxon | method | p/q).

#### Data to Improve Publication Impact

- **Bomb calorimetry of energy harvest:** Quantify energy extraction (kcal per gram) by ex vivo fecal microbiota under no-carbohydrate, RDC, and SDC conditions to test whether carbohydrate type and obesity status affect microbial energy harvest (aligns with Aim 2).
- **Energy harvest with physiologically relevant substrates:** Apply bomb calorimetry to ex vivo cultures supplied with pre-digested carbohydrates from a simulated gastrointestinal digestion model (e.g., SHIME or similar) to better approximate in vivo substrate availability and strengthen translational interpretation.
- **Mechanistic validation of *Fusicatenibacter*:** Measure SCFA production by *Fusicatenibacter* reference strains (e.g., *F. saccharivorans*, *F. faecihominis*) when grown on RDC vs SDC alone, and with acetate or Bifidobacterium feeder co-culture, to confirm cross-feeding and carbohydrate-utilization roles inferred from the integrated microbiome–SCFA analyses.

## Follow on Studies

- **Mechanistic follow-on:** Link ex vivo SCFA and taxa to in vivo outcomes (glycemia, weight, inflammation); validate *Fusicatenibacter* and other SDC-utilizers in gnotobiotic or cultured systems; test synbiotic (SDC + isolates) in preclinical or pilot human studies.
- **Larger/clinical:** Replicate in larger adolescent cohorts; test SDC interventions with microbiome and SCFA as mediators; define predictors of “SCFA responder” for trial stratification.

## Next Steps

### Near Term (3–6 months)

**1. Complete Aim 2 (Energy harvest quantification)**

- Perform bomb calorimetry on banked ex vivo culture samples (no-carb, RDC, SDC) from all subjects
- Analyze energy extraction patterns by obesity status and carbohydrate type
- Integrate energy harvest data with SCFA production and microbiome composition to test whether high-SCFA/low-energy phenotypes correlate with specific taxa
- Target completion: 3 months from sample retrieval and calorimetry setup

**2. Advance Aim 3 (SDC utilizer characterization) with taxa-focused approach**

- Obtain and culture reference strains:
  - *Fusicatenibacter saccharivorans* (DSM 26062) and *F. faecihominis* (DSM 113288)
  - *Bifidobacterium adolescentis* and *B. longum* reference strains
- Monoculture growth assays: Test SCFA production (acetate, lactate, butyrate, propionate) by each strain on RDC vs SDC substrates
- Co-culture experiments: *Bifidobacterium* + *Fusicatenibacter* to test acetate/lactate cross-feeding and enhanced butyrate production with SDC
- Compare reference strain SCFA profiles to ex vivo community patterns from responder vs non-responder subjects
- Target completion: 4–6 months (dependent on strain acquisition and culture optimization)

**3. Manuscript and presentation preparation**

- Draft full manuscript with current data (SCFA, microbiome, integrated analyses)
- Incorporate Aim 2 results as they become available; note Aim 3 progress in Discussion
- Prepare structured abstract and presentation for target meetings (DDW, PAS, Obesity Week, AGA, ASM Microbe)
- Submit abstract(s) for upcoming conference deadlines (typically 3–4 months in advance)
- Internal Abbott-Lurie presentation to review findings and discuss next-phase studies
- Target: Abstract submission within 6–8 weeks; manuscript draft for internal review within 3 months

### Mid-Term (6–18 months)

**1. Pilot dietary intervention: RDC vs SDC supplementation in adolescents**

*Rationale:* Translate ex vivo findings to in vivo setting; test whether SDC preferentially enriches *Fusicatenibacter*, *Bifidobacterium*, and other SCFA producers; identify microbiome/SCFA predictors of metabolic response.

**Key Questions:**

- Does SDC supplementation increase abundance of *Fusicatenibacter* and *Bifidobacterium* compared to RDC or control diet?
- Do baseline microbiome composition (presence/abundance of *Fusicatenibacter*, *Bifidobacterium*) and SCFA production capacity (ex vivo assay at baseline) predict glycemic, weight, or inflammatory responses to SDC?
- Are changes in *Fusicatenibacter* and *Bifidobacterium* abundance during intervention associated with responder vs non-responder phenotypes?

**Study Design:**

- Randomized, parallel-arm pilot study (n = 30–45 adolescents with overweight/obesity, ages 12–18)
- Three arms: (1) Control diet, (2) RDC supplementation (e.g., 30g/day maltodextrin or rapidly digestible starch), (3) SDC supplementation (e.g., 30g/day resistant starch type 2 or slowly digestible starch)
- Duration: 8–12 weeks intervention with 4-week follow-up
- Timepoints: Baseline, 4 weeks, 8 weeks (end of intervention), 12 weeks (follow-up)

**Primary Outcomes:**

- Change in fasting glucose and postprandial glucose (2-hour OGTT or continuous glucose monitoring)
- Change in HOMA-IR (fasting insulin and glucose)
- Change in body weight and BMI z-score

**Secondary Outcomes:**

- Fecal microbiome composition (16S and/or shotgun metagenomics): relative and absolute abundance of *Fusicatenibacter*, *Bifidobacterium*, other DA taxa from ex vivo study
- Fecal and plasma SCFA concentrations (acetate, propionate, butyrate)
- Inflammatory markers (hsCRP, IL-6, TNF-α)
- Dietary adherence and tolerability (food diaries, GI symptom questionnaires)
- Exploratory: breath hydrogen/methane, serum metabolomics, gut barrier markers (zonulin, LPS-binding protein)

**Methods and Measurements:**

- Recruitment: Lurie Children's Hospital adolescent weight management clinic
- Dietary intervention: Provided supplements (blinded sachets) mixed into food/beverages daily; monthly dietitian counseling to maintain isocaloric background diet
- Sample collection: Fasting blood, fecal sample (home collection kit, -80°C storage within 24h)
- Microbiome: 16S V4 sequencing + qPCR for absolute quantification; consider shotgun metagenomics for functional pathways
- SCFA: GC-MS PFBBr derivatization (same protocol as ex vivo study)
- Ex vivo fermentation assay at baseline: Use subject's baseline fecal sample in RDC/SDC culture to define "high" vs "low" SCFA producer phenotype for predictor analysis
- Statistical analysis: Mixed-effects models for repeated measures; correlation of baseline microbiome/ex vivo SCFA with outcome changes; responder definition (e.g., >10% reduction in HOMA-IR or glucose AUC)

**Timeline:** 12–18 months (3 months planning/IRB, 6–9 months enrollment/intervention, 3 months analysis/manuscript)

**2. Additional mid-term studies (concept stage)**

*Study of SDC to enhance GLP-1 agonist efficacy in adolescents*

- Rationale: GLP-1 agonists (semaglutide, liraglutide) are increasingly used for adolescent obesity; microbiome and SCFA modulation may enhance glycemic control, weight loss, or reduce GI side effects
- Design: Randomized add-on trial of SDC vs placebo in adolescents initiating GLP-1 agonist therapy
- Outcomes: Weight, glucose control, GI tolerability, microbiome/SCFA changes
- Sample size: n = 40–60 (pilot/feasibility)
- Timeline: 18–24 months (dependent on GLP-1 agonist use in clinic population)

*Study of SDC effects after transitioning off GLP-1 agonists*

- Rationale: Weight regain is common after GLP-1 agonist discontinuation; SDC may sustain microbiome/metabolic benefits and reduce rebound
- Design: Observational cohort or randomized maintenance trial (SDC vs control) in adolescents discontinuing GLP-1 agonist
- Outcomes: Weight trajectory, glucose, microbiome persistence of beneficial taxa, SCFA
- Timeline: 18–24 months (requires cohort of adolescents completing GLP-1 therapy)

*Study of SDC in adolescents with metabolic dysfunction-associated steatohepatitis (MASH)*

- Rationale: SCFA (especially butyrate, propionate) improve gut barrier, reduce endotoxemia and hepatic inflammation in preclinical NASH/MASH models; SDC may benefit pediatric MASH
- Design: Pilot RCT of SDC vs control in adolescents with biopsy-proven or MRI-confirmed MASH
- Outcomes: Hepatic fat (MRI-PDFF), ALT, inflammatory markers, microbiome/SCFA, histology (if repeat biopsy feasible)
- Collaboration: Pediatric hepatology, partnership with NASH clinical research network
- Timeline: 24+ months (requires specialty recruitment, imaging, potential biopsy)

### Long-Term Vision (2–5 years)

**1. Precision nutrition trials stratified by microbiome/SCFA phenotype**

- Multi-site RCTs of SDC ± synbiotic (*Fusicatenibacter*, *Bifidobacterium*, or consortia) in adolescents with obesity
- Enrollment stratified or enriched for baseline microbiome features (e.g., presence/absence of *Fusicatenibacter*, ex vivo high vs low SCFA producer)
- Adaptive trial designs: Early futility or efficacy signals based on microbiome/SCFA biomarkers
- Integration with digital health (CGM, activity trackers, microbiome monitoring)
- Outcomes: Weight, glycemia, cardiometabolic risk, quality of life, microbiome/SCFA durability

**2. Mechanistic and translational studies**

- Gnotobiotic mouse models: Colonize with defined consortia (*Fusicatenibacter* + *Bifidobacterium* + butyrate producers) ± SDC diet; measure weight, glucose, inflammation, intestinal barrier, hepatic metabolism
- Human organoid or epithelial co-culture models: Test SCFA (from ex vivo cultures or pure standards) effects on enterocyte glucose transport, GLP-1 secretion, barrier integrity
- Multi-omics integration: Metagenomics, metatranscriptomics, metabolomics (fecal, serum, urine) to map carbohydrate -> microbiome -> metabolite -> host pathways
- Collaboration with systems biology and computational groups for predictive modeling

**3. Extension to other populations and endpoints**

- Life stages: Adults with obesity/metabolic syndrome; younger children (school-age); pregnant individuals with gestational diabetes
- Cardiometabolic endpoints: Blood pressure, lipid profiles, arterial stiffness, cardiovascular events (long-term follow-up)
- Populations: Diverse race/ethnicity, socioeconomic backgrounds; international cohorts for microbiome variability
- Combination therapies: SDC + pharmacotherapy (metformin, SGLT2i, GLP-1 agonists); SDC + behavioral interventions

**4. Regulatory and commercialization pathway**

- If efficacy demonstrated: Develop medical food or dietary supplement product with standardized SDC composition
- Engage FDA (medical food vs drug classification) and payers (coverage for obesity/diabetes management)
- Partner with Abbott Nutrition for formulation, manufacturing, and distribution
- Companion diagnostic: Microbiome test to identify likely responders (precision nutrition stratification)

