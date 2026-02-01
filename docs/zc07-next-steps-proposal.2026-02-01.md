---
title: "Abbott-Lurie Children's Collaboration ZC07"
subtitle: "Review and Next Steps Proposal"
date: "February 2, 2026"
mainfont: "Roboto"
fontsize: 11pt
linestretch: 1
toc: true
geometry: "top=1in, bottom=1in, left=0.5in, right=0.5in"
output: pdf_document
header-includes:
  - \usepackage{fontspec}
  - \setmainfont[Scale=1.0]{Roboto}
  - \newfontfamily\fallbackfont{DejaVu Sans}
  - \usepackage{newunicodechar}
  - \newunicodechar{→}{{\fallbackfont →}}
  - \newunicodechar{α}{{\fallbackfont α}}
  - \newunicodechar{×}{{\fallbackfont ×}}
  - \defaultfontfeatures{Ligatures=TeX}
  - \usepackage{graphicx}
  - \usepackage{titlesec}
  - \titlespacing*{\subsubsection}{0pt}{3.25ex plus 1ex minus .2ex}{1.5ex plus .2ex}
  - \usepackage{wrapfig}
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhf{}
  - \fancyhead[R]{\textit{Abbott-Lurie Children's Collaboration ZC07, Updated 2026-02-02}}
  - \fancyfoot[C]{\thepage}
  - \fancypagestyle{plain}{\fancyhf{}\fancyfoot[C]{\thepage}\renewcommand{\headrulewidth}{0pt}}
---

*** 
## Introduction and Purpose 

This document provides a comprehensive review of the ZC07 Abbott-Lurie Children's Collaboration, summarizing progress to date and outlining a strategic roadmap for near-term, mid-term, and long-term activities to consider. The objectives of this documents are to:

1. **Review project goals and milestones:** Summarize the original aims and assess progress toward each objective.
2. **Identify completed and ongoing work:** Document analyses and experiments that have been finalized, and clarify remaining tasks.
3. **Address challenges:** Describe barriers encountered in completing outstanding aims and propose solutions.
4. **Outline dissemination strategy:** Present a structured abstract for conference presentation and a publication outline with key content, figures, and tables.
5. **Define mechanistic follow-on studies:** Identify experiments required to strengthen publication impact by elucidating mechanisms underlying observed phenotypes.
6. **Propose mid- and long-term investigations:** Articulate a research trajectory that extends current findings toward translational and clinical applications. 

***

## Project Background and Progress to Date 

### Significance

Obesity continues to increase as a public health emergency, with the origins of most adult obesity rooted in childhood. Effective clinical approaches are urgently needed to prevent or reverse childhood obesity. The gut-associated microbiome is an established central factor in energy harvest, hepatic function, insulin sensitivity, and adipose tissue homeostasis, making it a critical target for obesity intervention strategies. A growing body of evidence indicates that carbohydrate *quality*—particularly digestibility—matters more than quantity for metabolic outcomes. Resistant starch, slowly digestible starch (SDS), and high-molecular-weight β-glucans have been shown to improve glycemic control, reduce visceral adiposity, and modulate the gut microbiota in ways that favor insulin sensitivity and SCFA production. Meta-analyses and clinical trials support the therapeutic potential of these substrates for metabolic syndrome and obesity management, though optimal dosing, individual variability, and long-term sustainability remain active areas of research.

This project builds on foundational work demonstrating that slowly digestible carbohydrates (SDC) reduce glucose excursions in healthy, insulin-resistant, and type 2 diabetic individuals by inducing slow, prolonged glucose release and lower postprandial glycemic responses. In type 2 diabetic patients, SDC-rich diets have been shown to reduce glycemic variability parameters (e.g., standard deviation, MAGE) by 17–23%, with these parameters correlating with HbA1c. Abbott investigators have shown in rodent  models of obesity that nutrition with SDC, in comparison to fast digestible carbohydrate (FDC) sources, reverses obesity-associated phenotypes including elevated body mass, insulin resistance, and systemic inflammation. These preclinical findings provide a mechanistic rationale for testing SDC in human populations, but interpersonal differences in SDC utilization by the childhood-associated human microbiota may not be fully predicted by murine models.

Lurie investigators and colleagues have demonstrated interpersonal variation in short-chain fatty acid (SCFA) production by the human gut microbiome from adolescents with obesity in response to *ex vivo* prebiotic exposure. This variation suggests that complex carbohydrate utilization by the microbiota differs between individuals and may determine who responds to SDC and other nutritional approaches to obesity. Understanding the compositional and metabolic responses of the childhood-associated microbiota to slow versus fast digestible carbohydrates can therefore inform future obesity-treatment trials and precision approaches to therapy—including the identification of likely responders and the development of synbiotic or dietary strategies tailored to microbiome phenotype. The ZC07 collaboration between Abbott Nutrition and Lurie Children's Hospital was established to measure this variation systematically and to identify childhood-associated organisms with facile utilization of SDC, with the long-term goal of enabling precision nutrition interventions for pediatric obesity.

### Project Goals

- Measure variation in responses of human gut-associated microbiome communities to fast (FDC/RDC) and slow digestible (SDC) carbohydrates.
- Identify childhood-associated organisms with facile utilization of SDC.
- Test interindividual variation of short-chain fatty acid (SCFA) production among fecal microbiota samples to slow and fast digestible carbohydrates.
- Measure energy harvest differences between human obesity-associated fecal microbiota (Aim 2).
- Isolate SDC bacterial utilizers using single-cell isolation toward synbiotic combinations (Aim 3).

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{./images/approach-overview.png}
\caption{\small Approach overview}
\end{figure}
### Analysis Progress

#### Microbiome Analysis
*Reference: microbiome_analysis/analysis.v4_zymo.html*

16S rRNA gene sequencing (Zymo pipeline) with phyloseq-based analysis; alpha diversity (Observed, Shannon, Simpson) at 0H and 48H; beta diversity (Bray-Curtis PCoA, PERMANOVA with subject strata); baseline clustering and trajectories; genus-level composition changes 0H->48H by carbohydrate type; expanded differential abundance (corncob, ANCOM-BC2, LEfSe, MaAsLin2, etc.) with 48H focus; FDR control throughout.

#### SCFA Metabolomics Analysis
*Reference: scfa_metabolomics/scfa-project2-analysis-improved-v2.html*

Targeted SCFA metabolomics (GC-MS/PFBBr: acetate, butyrate, propionate, 5-aminovalerate, succinate); summary statistics by group and carbohydrate type; SCFA concentrations over time; delta (48h-0h) response magnitude; two-way ANOVA (group × carbohydrate) on deltas; carbohydrate type effects (combined and stratified by Control/Case); mixed-effects models for SCFA trajectories. Key outcomes: carbohydrate type and time drive SCFA production; RDC and SDC both increase SCFAs vs no-carbohydrate; limited Case vs Control differences; pronounced interpersonal variation in SCFA response.

#### Integrated Analysis
*Reference: integrated/integrated-scfa-microbiome.html*

Single workflow combining SCFA and microbiome with shared metadata; SCFA delta and mixed-effects results; alpha/beta diversity; composition and differential abundance; taxon–SCFA associations (e.g., *Fusicatenibacter* vs butyrate/propionate); responder vs non-responder (median-delta split for butyrate/propionate) and composition comparisons; absolute abundance (qPCR-based) for focal taxa; correlation and regression for taxa–SCFA relationships.

### Progress to Date

This project has successfully completed comprehensive *ex vivo* analysis of fecal microbiota from adolescents with and without obesity, examining metabolic and compositional responses to carbohydrate substrates of varying digestibility. The integrated analysis combines targeted metabolomics, 16S rRNA gene sequencing, and advanced statistical modeling to provide mechanistic insights into interpersonal variation in carbohydrate fermentation capacity.

#### Major Findings

\begin{wrapfigure}{r}{0.6\textwidth}
\vspace{-2\baselineskip}
\centering
\includegraphics[width=\linewidth]{./images/key-findings-summary.png}
\vspace{-25pt}
\caption{\small Key Findings To Date}
\vspace{0\baselineskip}
\end{wrapfigure}

1. **Carbohydrate type and time are primary determinants of SCFA production.** Both rapid digestible carbohydrate (RDC) and slow digestible carbohydrate (SDC) elicited significant increases in butyrate, propionate, and acetate concentrations over 48 hours compared to no-carbohydrate controls. Mixed-effects models accounting for repeated measures confirmed robust main effects and carbohydrate × timepoint interactions, demonstrating that carbohydrate structure and digestibility shape microbial metabolic output.

\begin{wrapfigure}{r}{0.6\textwidth}
\vspace{-2\baselineskip}
\centering
\vspace{0\baselineskip}
\end{wrapfigure}

2. **Obesity status does not substantially impair *ex vivo* fermentation capacity.** While minor baseline differences in SCFA levels were observed, the overall magnitude and trajectory of SCFA responses to RDC and SDC were comparable between adolescents with obesity (Case) and healthy-weight controls. Limited significant group main effects or interactions in mixed-effects models indicate that the fundamental microbial capacity to ferment these substrates is preserved in the context of obesity in this cohort.

3. **Pronounced interpersonal variation in SCFA production supports precision nutrition approaches.** Subject-to-subject variability in SCFA response magnitude was substantial, with some individuals consistently classified as high or low producers across carbohydrate conditions. This heterogeneity persisted even after accounting for obesity status and carbohydrate type, underscoring that individual microbiome composition—rather than obesity phenotype alone—drives metabolic responsiveness. Responder vs non-responder stratification (median-delta split for butyrate and propionate) enabled targeted characterization of microbiome features associated with differential SCFA production.

4. **Specific taxa associate with SCFA concentrations and responder phenotypes.** Differential abundance analyses (ANCOM-BC2, MaAsLin2, LEfSe) identified multiple taxa that differ significantly between carbohydrate conditions, timepoints, and responder groups. Notably, *Fusicatenibacter* (Lachnospiraceae) exhibited positive associations with butyrate and propionate concentrations, particularly under SDC conditions at 48 hours. Correlation and regression analyses confirmed monotonic relationships between *Fusicatenibacter* relative and absolute abundance (qPCR-derived) and SCFA levels, implicating this genus as a candidate SDC utilizer and butyrate producer. *Bifidobacterium* abundance increased with both RDC and SDC (more so with SDC), consistent with its known role as a primary degrader of complex carbohydrates and cross-feeder for secondary fermenters.

5. **Microbiome composition shifts in response to carbohydrate challenge.** Alpha diversity (Observed richness, Shannon, Simpson indices) showed modest changes between baseline (0H) and 48H, with some carbohydrate-specific effects. Beta diversity analyses (Bray-Curtis PCoA, PERMANOVA) revealed significant effects of carbohydrate type, timepoint, and their interactions on community structure, with subject ID as a key stratification variable capturing baseline compositional differences. Baseline clustering and trajectory analyses demonstrated that individuals exhibit distinct community shifts in response to RDC vs SDC, further supporting personalized metabolic responses.

6. **Integration of SCFA and microbiome data identifies taxon–metabolite linkages.** By incorporating SCFA concentrations as continuous predictors in differential abundance models, taxa positively or negatively associated with butyrate and propionate were identified. These taxon–SCFA associations provide mechanistic hypotheses for how specific community members contribute to overall fermentation output. The responder stratification revealed compositional differences: high SCFA responders had enriched abundance of putative butyrate producers (e.g., *Fusicatenibacter*, *Roseburia*) and primary fermenters (e.g., *Bifidobacterium*), while low responders showed reduced representation of these taxa.

#### Implications
These findings establish that carbohydrate quality, rather than obesity status per se, is a modifiable determinant of microbial SCFA production in adolescents. The substantial interpersonal variation and identification of candidate taxa (*Fusicatenibacter*, *Bifidobacterium*) that associate with SCFA responsiveness provide a foundation for precision nutrition strategies. Future interventions may be optimized by: (1) selecting individuals with baseline microbiome profiles enriched in SDC-utilizing taxa, (2) co-administering targeted probiotics or synbiotics alongside SDC to enhance fermentation, or (3) using *ex vivo* SCFA assays as companion diagnostics to predict *in vivo* responsiveness to dietary carbohydrate interventions.

#### Gaps in Translation
Several gaps limit translation of these *ex vivo* findings to human health outcomes. First, elevations in fecal SCFA do not necessarily correlate with plasma SCFA; systemic availability and metabolic fate of microbially produced SCFAs remain to be established in this population. Second, microbiome changes observed under specific carbohydrate conditions *ex vivo* may not reflect what occurs *in vivo*, where substrate availability, transit time, and host factors differ. Third, elevated SCFA production—whether fecal or plasma—may or may not associate with meaningful clinical outcomes such as improved insulin sensitivity, HbA1c, or weight loss; these relationships have not been tested in this cohort. Fourth, *Fusicatenibacter* and *Bifidobacterium*, though associated with SCFA production *ex vivo*, may or may not correlate with these human outcomes, and whether their presence or abundance in the baseline microbiota predicts dietary response to SDC is unknown.

Additional gaps to consider include: **gastrointestinal transit and digestion**—the carbohydrates in this study were not exposed to different segments of the gastrointestinal tract (stomach, small intestine) where they may be partially digested or modified before reaching the colon; *ex vivo* fermentation therefore uses substrates that may differ from those available to the colonic microbiota *in vivo*; **dose, formulation, and duration**—*ex vivo* conditions use fixed substrates and short incubations (48–72 h), whereas *in vivo* dose, SDC type (e.g., resistant starch vs slowly digestible starch), food matrix, and intervention length (weeks to months) will affect substrate delivery and microbiome adaptation; **dietary context**—*ex vivo* tests single substrates, but *in vivo* the background diet (fiber, fat, other carbohydrates) and whether SDC is substituted or added will influence response; **host factors**—medication use (e.g., metformin, GLP-1 agonists), genetics, gut motility, and intestinal transit are not captured *ex vivo* but can modify *in vivo* microbiome and metabolic outcomes; **generalizability**—this cohort is a single adolescent population from one site, and results may not extend to other ages, ethnicities, or obesity severity; and **causality**—*ex vivo* taxon–SCFA associations are correlative, and causal roles in health outcomes require intervention or mechanistic studies. Closing these gaps requires a direct human study (e.g., a pilot dietary intervention) in which SDC or control supplementation is administered, microbiome and SCFA are measured over time, and glycemic, metabolic, and weight outcomes are assessed to test whether *ex vivo* phenotypes predict *in vivo* response.

***

## Steps in Progress

\begin{figure}[h]
\centering
\includegraphics[width=0.5\textwidth]{./images/cross-feeding-diagram.png}
\vspace{-10pt}
\caption{\small Cross-feeding mechanism}
\end{figure}

- **Aim 2 (energy harvest):** *Ex vivo* exposure of fecal samples to SDC and FDC with timepoint collection and bomb calorimetry (kcal per gram) not yet completed.
- **Aim 3 (SDC utilizers):** Single-cell capture, arraying, and growth assays for SDC vs FDC utilizers from high-SCFA / low–energy-extraction donors not yet completed.
    - **Reference strain acquisition (*Fusicatenibacter*):** Fully sequenced type strains are being procured from the DSMZ culture collection. Two strains have been identified: *F. saccharivorans* (DSM 26062; type strain HT03-11), originally isolated from human feces, with documented saccharolytic capacity for glucose, arabinose, xylose, cellobiose, and starch, and established relevance to carbohydrate metabolism and SCFA production within Lachnospiraceae; and *F. faecihominis* (DSM 113288), also of human fecal origin, with a complete genome sequence available for comparative genomic and functional validation studies.
    - **Reference strain acquisition (*Bifidobacterium*):** Given the observed enrichment of *Bifidobacterium* under both RDC and SDC conditions (with greater increases under SDC), procurement of representative strains is warranted. *B. adolescentis* and *B. longum* are established primary fermenters and cross-feeders for butyrate-producing taxa: these species ferment oligosaccharides and starch to lactate and acetate, which are subsequently converted to butyrate by lactate-utilizing organisms (e.g., *Anaerostipes*, *Eubacterium hallii*, *Roseburia*, *Faecalibacterium prausnitzii*). Additionally, *Bifidobacterium* species release partially hydrolyzed oligosaccharides that serve as substrates for select butyrate producers. Reference strains of these species will enable mechanistic and synbiotic co-culture experiments in conjunction with *Fusicatenibacter*.
- **Analysis:** Finalizing responder/non-responder definitions and validation; any remaining cross-dataset QC and narrative for publication.

**Challenges in completing remaining goals:** Resource and timeline constraints for bomb calorimetry and single-cell isolation workflows; need to align sample availability and lab capacity for Aims 2–3; defining responder thresholds (median-delta vs clinical anchors) for publication; coordinating integrated analyses across SCFA and microbiome pipelines.

***

## Near Term Work 

Near-term priorities center on disseminating current findings and strengthening the manuscript. Below are a draft conference abstract, a publication outline with planned content and figures, and a list of additional data that would improve publication impact.

### Abstract (Sample)
 *Microbiome Composition and Short-Chain Fatty Acid Response to Slow and Rapid Digestible Carbohydrates: An Ex Vivo Study in Adolescents With and Without Obesity*

**Introduction.** Childhood obesity is a major driver of adult obesity, and the gut microbiome is a central mediator of metabolic health. Carbohydrate quality—particularly digestibility—may matter more than quantity for metabolic outcomes. Interpersonal variation in microbial short-chain fatty acid (SCFA) production may determine who benefits from slowly digestible carbohydrates (SDC). We aimed to characterize *ex vivo* SCFA production and microbiome composition in adolescents with and without obesity in response to different carbohydrate substrates.

**Methods.** Fecal samples from adolescents with obesity (Case) and healthy-weight controls (Control) were incubated *ex vivo* under no-carbohydrate, rapid digestible (RDC), and slow digestible (SDC) conditions. SCFAs (acetate, butyrate, propionate, and related metabolites) were quantified by gas chromatography–mass spectrometry after PFBBr derivatization; microbiome composition was assessed by 16S rRNA gene sequencing. Analyses included linear mixed-effects models for SCFA trajectories, alpha and beta diversity with PERMANOVA, differential abundance (ANCOM-BC2, MaAsLin2), taxon–SCFA associations, and responder vs non-responder stratification (median-delta SCFA).

**Results.** Carbohydrate type and incubation time were the primary drivers of SCFA production; both RDC and SDC increased SCFAs versus no-carbohydrate controls with SDC driving the highest levels of butyrate. Obesity status did not affect fermentation capacity. Pronounced interpersonal variation in SCFA response was observed. Specific taxa (e.g., *Fusicatenibacter*) were associated with SCFA concentrations; responder stratification enabled targeted characterization of microbiome composition differences.

**Conclusions.** *Ex vivo* fermentation capacity for RDC and SDC is preserved in adolescents with obesity, with substantial interindividual variation in SCFA production. Taxon–SCFA associations and responder phenotypes support precision nutrition strategies and future trials stratified by microbiome/SCFA phenotype.

**Target audience/meetings:** Pediatric obesity, nutrition, microbiome (e.g., DDW, PAS, Obesity Week, AGA, ASM Microbe); Abbott-internal and Lurie–Abbott collaboration meetings.

### Publication 

- **Title options (examples):**
  1. → *Microbiome Composition and Short-Chain Fatty Acid Response to Slow and Rapid Digestible Carbohydrates: An Ex Vivo Study in Adolescents With and Without Obesity*
  2. *Interpersonal Variation in Short-Chain Fatty Acid Production by Adolescent Gut Microbiota in Response to Digestible Carbohydrate Substrates*
  3. *Carbohydrate Digestibility Drives Ex Vivo SCFA Production Independent of Obesity Status in Adolescent Fecal Microbiota*
  4. *Fusicatenibacter and Bifidobacterium Associate with SCFA Production Capacity in Response to Digestible Carbohydrates: Implications for Precision Nutrition in Pediatric Obesity*
  5. *Ex Vivo Fermentation Capacity of Adolescent Gut Microbiota: Carbohydrate Type Supersedes Obesity Status as a Determinant of SCFA Production*
- **Structured abstract format:** Background, Methods, Results, Conclusions (see sample abstract above).
- **Introduction:** Childhood obesity, microbiome as metabolic mediator, SDC/FDC rationale, Lurie–Abbott aims, interpersonal variation and precision nutrition.
- **Methods:** Subject/sample description; *ex vivo* culture conditions; SCFA (extraction, derivatization, GC-MS); 16S sequencing and processing; alpha/beta diversity, differential abundance, mixed-effects and taxon–SCFA analyses; responder definition; software and FDR.
- **Results (content we have):** Sample summary; SCFA by group and carbohydrate type; SCFA over time and delta analyses; mixed-effects results; alpha/beta diversity; composition and differential abundance; taxon–SCFA (e.g., *Fusicatenibacter*) and responder vs non-responder composition.
- **Results (to fill in):** Bomb calorimetry (if Aim 2 completed); isolate growth data (if Aim 3 used); any additional validation or sensitivity analyses.
- **Discussion:** Interpretation of carbohydrate and time effects; interpersonal variation and precision nutrition; taxa–SCFA and responder findings; limitations (*ex vivo*, single cohort, timepoints); future directions (Aim 2/3, mechanisms, trials).
- **Major figures/tables planned:** (1) Study design/sampling schema; (2) SCFA by carbohydrate type and time (e.g., line + delta boxplots); (3) Alpha/beta diversity by group, carbohydrate, timepoint; (4) Composition (e.g., genus bar/alluvial) 0H vs 48H; (5) Differential abundance (forest or volcano); (6) Taxon–SCFA correlations (e.g., *Fusicatenibacter* vs butyrate/propionate); (7) Responder vs non-responder composition; (8) Summary tables: mixed-effects SCFA, PERMANOVA, DA summary (taxon | method | p/q).

### Data to Improve Publication Impact

- **Bomb calorimetry of energy harvest:** Quantify energy extraction (kcal per gram) by *ex vivo* fecal microbiota under no-carbohydrate, RDC, and SDC conditions to test whether carbohydrate type and obesity status affect microbial energy harvest (aligns with Aim 2).
- **Energy harvest with physiologically relevant substrates:** Apply bomb calorimetry to *ex vivo* cultures supplied with pre-digested carbohydrates from a simulated gastrointestinal digestion model (e.g., SHIME or similar) to better approximate *in vivo* substrate availability and strengthen translational interpretation.
- **Mechanistic validation of *Fusicatenibacter*:** Measure SCFA production by *Fusicatenibacter* reference strains (e.g., *F. saccharivorans*, *F. faecihominis*) when grown on RDC vs SDC alone, and with acetate or Bifidobacterium feeder co-culture, to confirm cross-feeding and carbohydrate-utilization roles inferred from the integrated microbiome–SCFA analyses.

***

## Follow on Studies for Consideration

- **Mechanistic follow-on:** Link *ex vivo* SCFA and taxa to *in vivo* outcomes (glycemia, weight, inflammation); validate *Fusicatenibacter* and other SDC-utilizers in gnotobiotic or cultured systems; test synbiotic (SDC + isolates) in preclinical or pilot human studies.
- **Larger/clinical:** Replicate in larger adolescent cohorts; test SDC interventions with microbiome and SCFA as mediators; define predictors of “SCFA responder” for trial stratification.

***

## Next Steps

### Near Term (3–6 months)

**1. Complete Aim 2 (Energy harvest quantification)**

- Perform bomb calorimetry on banked *ex vivo* culture samples (no-carb, RDC, SDC) from all subjects
- Analyze energy extraction patterns by obesity status and carbohydrate type
- Integrate energy harvest data with SCFA production and microbiome composition to test whether high-SCFA/low-energy phenotypes correlate with specific taxa
- Target completion: 3 months from sample retrieval and calorimetry setup


\begin{figure}[h]
\centering
\includegraphics[width=0.7\linewidth]{./images/near-mid-term-goals.png}
\vspace{-10pt}
\caption{\small Near- and mid-term goals}
\end{figure}

**2. Advance Aim 3 (SDC utilizer characterization) with taxa-focused approach**

- Obtain and culture reference strains:
  - *Fusicatenibacter saccharivorans* (DSM 26062) and *F. faecihominis* (DSM 113288)
  - *Bifidobacterium adolescentis* and *B. longum* reference strains
- Monoculture growth assays: Test SCFA production (acetate, lactate, butyrate, propionate) by each strain on RDC vs SDC substrates
- Co-culture experiments: *Bifidobacterium* + *Fusicatenibacter* to test acetate/lactate cross-feeding and enhanced butyrate production with SDC
- Compare reference strain SCFA profiles to *ex vivo* community patterns from responder vs non-responder subjects
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

*Rationale:* Translate *ex vivo* findings to *in vivo* setting; test whether SDC preferentially enriches *Fusicatenibacter*, *Bifidobacterium*, and other SCFA producers; identify microbiome/SCFA predictors of metabolic response.

**Key Questions:**

- Does SDC supplementation increase abundance of *Fusicatenibacter* and *Bifidobacterium* compared to RDC or control diet?
- Do baseline microbiome composition (presence/abundance of *Fusicatenibacter*, *Bifidobacterium*) and SCFA production capacity (*ex vivo* assay at baseline) predict glycemic, weight, or inflammatory responses to SDC?
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

- Fecal microbiome composition (16S and/or shotgun metagenomics): relative and absolute abundance of *Fusicatenibacter*, *Bifidobacterium*, other DA taxa from *ex vivo* study
- Fecal and plasma SCFA concentrations (acetate, propionate, butyrate)
- Inflammatory markers (hsCRP, IL-6, TNF-α)
- Dietary adherence and tolerability (food diaries, GI symptom questionnaires)
- Exploratory: breath hydrogen/methane, serum metabolomics, gut barrier markers (zonulin, LPS-binding protein)

**Methods and Measurements:**

- Recruitment: Lurie Children's Hospital adolescent weight management clinic
- Dietary intervention: Provided supplements (blinded sachets) mixed into food/beverages daily; monthly dietitian counseling to maintain isocaloric background diet
- Sample collection: Fasting blood, fecal sample (home collection kit, -80°C storage within 24h)
- Microbiome: 16S V4 sequencing + qPCR for absolute quantification; consider shotgun metagenomics for functional pathways
- SCFA: GC-MS PFBBr derivatization (same protocol as *ex vivo* study)
- *Ex vivo* fermentation assay at baseline: Use subject's baseline fecal sample in RDC/SDC culture to define "high" vs "low" SCFA producer phenotype for predictor analysis
- Statistical analysis: Mixed-effects models for repeated measures; correlation of baseline microbiome/*ex vivo* SCFA with outcome changes; responder definition (e.g., >10% reduction in HOMA-IR or glucose AUC)

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
- Enrollment stratified or enriched for baseline microbiome features (e.g., presence/absence of *Fusicatenibacter*, *ex vivo* high vs low SCFA producer)
- Adaptive trial designs: Early futility or efficacy signals based on microbiome/SCFA biomarkers
- Integration with digital health (CGM, activity trackers, microbiome monitoring)
- Outcomes: Weight, glycemia, cardiometabolic risk, quality of life, microbiome/SCFA durability

**2. Mechanistic and translational studies**

- Gnotobiotic mouse models: Colonize with defined consortia (*Fusicatenibacter* + *Bifidobacterium* + butyrate producers) ± SDC diet; measure weight, glucose, inflammation, intestinal barrier, hepatic metabolism
- Human organoid or epithelial co-culture models: Test SCFA (from *ex vivo* cultures or pure standards) effects on enterocyte glucose transport, GLP-1 secretion, barrier integrity
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

