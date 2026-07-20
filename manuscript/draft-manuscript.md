# DRAFT MANUSCRIPT — v0.1 (initial rough draft)

**Project:** Abbott–Lurie ZC07, Project 2
**Generated:** 2026-07-19
**Status:** ROUGH FIRST DRAFT for author revision. Not submission-ready.

---

## Title

**Incubation Time, Carbohydrate Condition, and Donor Identity Shape Ex Vivo Short-Chain Fatty Acid Accumulation by Adolescent Fecal Microbiota**

*Running title:* Ex vivo SCFA responses of adolescent fecal microbiota

**Authors:** ⟦PLACEHOLDER: author list, affiliations, ORCIDs, corresponding author⟧

**Conflict of interest statement:** ⟦PLACEHOLDER: Abbott funding/affiliation disclosure. Note that the originating rodent study (Plaza-Díaz 2022) and the related pediatric crossover trial (Gillen 2021, NCT03185884) are Abbott-affiliated; disclose transparently.⟧

---

## Abstract

**Background.** Carbohydrate structure may influence gut microbial fermentation, but fermentation responses to the same substrate vary between individuals. Whether obesity status or donor-level microbial features explain this variation in adolescents is uncertain.

**Methods.** Fecal microbiota from 16 adolescent donors (8 healthy-weight, 8 with obesity) were incubated ex vivo under three paired conditions — no added carbohydrate, a rapid-digesting carbohydrate formulation (RDC), and a slow-digesting carbohydrate formulation (SDC) — and sampled at 0 and 48 h. Acetate, propionate, butyrate, succinate and 5-aminovalerate were quantified by pentafluorobenzyl bromide (PFBBr) derivatization with negative-chemical-ionization GC-MS. Communities were profiled by 16S rRNA V3–V4 sequencing and scaled using total bacterial 16S qPCR. Analyses respected the donor → A/B aliquot → culture-well hierarchy, with donor as the independent biological unit.

**Results.** In donor-aware 48-h endpoint models, carbohydrate condition was associated with acetate (F(2,28)=16.11, q=6.63×10⁻⁵), propionate (F(2,28)=25.16, q=1.67×10⁻⁶), butyrate (F(2,28)=8.58, q=0.00371) and succinate (F(2,28)=5.46, q=0.0298), but not 5-aminovalerate (q=0.784); no obesity-group or group×carbohydrate term survived adjustment. Acetate, propionate and butyrate accumulated over 48 h in every condition, including no added carbohydrate. Adjusted pairwise change contrasts were analyte-specific and did not favor SDC: for the prespecified butyrate contrast, mean 0–48 h change was 7.95 mM under SDC versus 6.92 mM under RDC (adjusted difference 1.03 mM, 95% CI −0.85 to 2.91, p=0.396). None of six obesity-group difference-in-change contrasts was significant; equivalence was not evaluated. Donor identity accounted for substantial repeatable variation (ICC 0.63–0.84; RDC–SDC Spearman ρ 0.83–0.96). In donor-aggregated PERMANOVA, incubation time (R²=0.239, p=0.0001) and obesity group (R²=0.026, p=0.0037) were associated with Bray–Curtis structure, whereas carbohydrate (R²=0.014, p=0.665) and group×carbohydrate (p=1.000) were not. An exploratory ANCOM-BC2 analysis identified higher relative *Fusicatenibacter* abundance under SDC than RDC at 48 h (q<0.1, exploratory); qPCR-scaled models did not establish a significant SDC-associated expansion, and the donor-level association between SDC-attributable *Fusicatenibacter* and butyrate change was positive but nonsignificant (ρ=0.477, permutation p=0.088, bootstrap 95% CI −0.239 to 0.892, N=14).

**Conclusions.** Carbohydrate condition was associated with endpoint concentrations of four fermentation metabolites at 48 h, but adjusted changes did not establish greater SCFA accumulation under SDC than RDC. Ex vivo accumulation varied by incubation time and, substantially, by donor. Obesity status did not significantly modify SCFA changes, although equivalence was not established. Community and taxon findings are hypothesis-generating.

**Keywords:** gut microbiota; short-chain fatty acids; ex vivo fermentation; adolescent obesity; slowly digestible carbohydrate; interpersonal variation

---

## 1. Introduction

Adolescent obesity is common and carries substantial cardiometabolic comorbidity [NCDRisC2024; GBD2025; Marcus2022]. A feature of this population that complicates intervention design is the breadth of *between-person* variation in response. Across adolescent cohorts spanning three intervention types, individual BMI change has ranged from −50.2% to +12.9%; within the lifestyle-intervention arm alone the range was −25.4% to +5.0%, while bariatric surgery accounted for the largest reductions [Ryder2019]. Comparable heterogeneity is seen in dietary and glycemic responses in adults [Zeevi2015; PREDICT1]. This variation is usually framed as a property of the host. It may also, in part, be a property of the resident microbial community, which differs between individuals both in composition and in how it metabolizes a given substrate. These are distinct phenotypes and are not interchangeable: a host metabolic response and a microbial community's fermentation output can diverge, and neither predicts the other by default. Whether the gut microbiome is a cause or a correlate of pediatric obesity remains contested [Sze2016; Meijnikman2018; Ridaura2013], and the present study is not designed to resolve that question.

Carbohydrates that escape host digestion are fermented by colonic bacteria to short-chain fatty acids (SCFAs), principally acetate, propionate and butyrate, alongside intermediates such as succinate and lactate [Louis2009; Koh2016]. Butyrate is largely a secondary product, formed in part by cross-feeding: several butyrogenic taxa consume acetate while producing butyrate [Duncan2002; Duncan2004; Belenguer2006]. A measured concentration therefore reflects **net accumulation** — production minus consumption and conversion — rather than metabolic flux [Sakata2019]. This distinction is not merely semantic. In humans, only about 36%, 9% and 2% of colonically delivered acetate, propionate and butyrate reach the systemic circulation [Boets2017], and fecal concentrations do not track modelled production rates whereas plasma concentrations do [Kirschner2025]. Accordingly, culture SCFA concentrations index microbial fermentation potential under defined exposure; they do not index colonic delivery, host exposure, or clinical benefit. Higher fecal SCFA is not straightforwardly favorable — in one large cohort it associated with lower diversity, greater intestinal permeability, obesity and hypertension [deLaCuestaZuluaga2018], a finding most easily reconciled if fecal SCFA partly reports absorption efficiency rather than production.

Carbohydrate digestibility class — rapidly digestible, slowly digestible, and resistant starch — is defined by the *rate of host enzymic hydrolysis*, with only the resistant fraction offering a guide to colonic delivery [Englyst1992]. The motivating observation for the present work comes from a rodent study in which a rapidly digestible carbohydrate source within a high-fat diet was replaced by a multi-component slower-digesting blend, with accompanying microbiome shifts and higher fecal acetate and propionate [PlazaDiaz2022]. Two features of that study bound what it can motivate. Its butyrate contrast was not significant, so it does not establish a butyrate expectation. And because the substituted diet also exchanged cellulose for resistant maltodextrin plus inulin/fructo-oligosaccharides at matched declared fiber, digestion rate and fermentable-substrate identity were varied together, so its SCFA findings cannot be attributed to digestion rate alone. Related human and in vitro work indicates that carbohydrate structure more often redistributes SCFA *among* analytes than raises them uniformly [Giuberti2013; Deehan2020], and that "more substrate" does not reliably yield "more SCFA" once background medium, pH, competition and cross-feeding are considered [Poppe2023; VanWehle2024].

An ex vivo fermentation assay applies substrate directly to a fecal community, deliberately removing host digestion, absorption, transit, mucosal surfaces and immune signaling [PerezBurillo2021; Isenring2023]. This isolates microbial response potential, and it also changes the independent variable: a formulation's host-digestibility class no longer governs what the community sees. That is not a subtle gap. Isomaltulose, which constitutes 26.4% of the slower-digesting blend used in the originating rodent work, showed small-intestinal digestibility of 95.5% and 98.8% and absorption of 93.6% and 96.1% in human ileostomy studies [Holub2010]; in vivo it therefore functions largely as a glycemic manipulation that reaches the colon minimally, whereas an ex vivo assay presents it intact. Static digestion protocols cannot repair this, since they are explicitly unsuitable for simulating digestion *kinetics* [Brodkorb2019]. Ex vivo systems remain well suited to a different question — how much donor communities differ in their intrinsic fermentation response to a defined substrate — and donor identity is a large source of variance in such assays [Reichardt2018; Holmes2020]. Evidence in adolescents specifically is sparse: a targeted search identified essentially no adolescent-specific ex vivo carbohydrate fermentation literature beyond a single donor-resolved study in adolescents with obesity, which lacked a healthy-weight comparator [Holmes2020].

We characterized donor-specific SCFA and microbial community responses of adolescent fecal microbiota to rapid- and slow-digesting carbohydrate formulations ex vivo, and evaluated whether these responses differed by obesity status or were associated with microbial community features.

Primary hypotheses were: (1) SCFA trajectories would vary across no-added-carbohydrate, RDC and SDC conditions; (2) SDC would produce a greater 0-to-48-h butyrate increase than RDC, with estimand (SDC 48 h − SDC 0 h) − (RDC 48 h − RDC 0 h); (3) obesity status would modify carbohydrate-specific SCFA change; and (4) SCFA responses would vary among donors.

Exploratory objectives were to test community-level associations with time, carbohydrate and obesity group; to identify carbohydrate-associated genera using relative and qPCR-scaled frameworks; to evaluate continuous donor-level taxon–SCFA relationships; and to describe median-defined response strata without treating them as validated phenotypes.

## 2. Materials and methods

**Study population and ethics**

[Original POMMS study and new JCI paper]. Samples were from individuals ages X to Y years. Donors were categorized as obese weight or healthy weight based on ____.  All subjects were not on therapy at the time of the sample collections used in this study.  Samples were deidentified. 

**Fecal collection and aliquoting**
Details of fecal sample collection and initial storage are provided in [references]. Samples from the original POMMS biorepository were shipped and stored frozen at -80 C until use in the studies herein. At the time of performing assays, samples were thawed at room temperature under anaerobic conditions. The samples were manually homogenized and split into samples required for biological and technical replicates with the remaining samples split into aliquots for follow up assays and refrozen at -80 C.

***Ex vivo* fermentation**

Biological replicate fecal inocula were incubated in technical duplicate anaerobically in pH balanced MiPro medium [reference] without peptone under three paired conditions: no added carbohydrate, RDC, and SDC. Cultures were sampled at 0 and 48 h. Within each aliquot and condition, culture wells were designated R1/R2 and S1/S2; these are culture-well replicates nested within aliquot × condition, not independent samples.

⟦PLACEHOLDER: exact RDC and SDC composition, source, and lot; carbohydrate dose (% w/v); MIPRO medium formulation and starting pH; inoculum concentration (% w/v); anaerobic atmosphere composition and method; incubation temperature; vessel type and working volume; agitation; plate allocation and randomization; termination and sample-storage method; background fermentable carbohydrate content of the medium.⟧

**SCFA quantification**

Acetate, propionate, butyrate, and succinate were quantified by the University of Chicago Duchossois Family Institute Host-Microbe Metabolomics Facility (DFI-HMMF) using GC-MS after pentafluorobenzyl bromide derivatization, as previously described [Haak2018], with facility-standard QC. We report SCFA **concentration** and **net accumulation** (0-to-48-h change). We avoid "production," because substrate disappearance and metabolic flux were not measured [Sakata2019].

**16S rRNA sequencing and processing**

Bacterial 16S rRNA gene community profiling was performed with the ZymoBIOMICS® Targeted Sequencing Service (Zymo Research, Irvine, CA; project zr24558; run dated 2025-08-13). DNA was extracted using the ZymoBIOMICS®-96 MagBead DNA Kit. The V3–V4 region was amplified and libraries prepared with the *Quick*-16S™ NGS Library Prep Kit and Zymo Research custom V3–V4 primers under real-time PCR cycle control to limit chimera formation; libraries were then equimolarly pooled. Provider positive controls (ZymoBIOMICS® Microbial Community Standard for extraction; Microbial Community DNA Standard for library preparation) and blank extraction and library-preparation negative controls were included. The pooled library was sequenced on an Illumina® NextSeq™ with a P1 reagent kit (600 cycles) and 30% PhiX spike-in. Unique amplicon sequence variants (ASVs) were inferred by the provider with DADA2, including removal of putative sequencing errors and chimeras [Callahan2016]; the delivered ASV table comprised 2,620 ASVs (mean length 406.7 bp).

Downstream handling used phyloseq [McMurdie2013], with relative abundances computed as within-sample proportions and analyses conducted at genus level.

**Total bacterial qPCR and qPCR-scaled genus abundance**

Total bacterial 16S load was quantified by Zymo Research as part of the same Targeted Sequencing Service using quantitative real-time PCR with the same V3–V4 primers used for library preparation. A plasmid standard containing one copy of the 16S gene (and one copy of the fungal ITS2 region) was prepared in 10-fold serial dilutions; gene copies per microliter of DNA extract were calculated from the standard curve and the PCR input volume. Provider-reported genome copies per microliter were derived by dividing gene copies by an assumed four 16S copies per genome; claim-facing analyses used the measured gene-copy values (`gene_copies_per_ul`), not the genome-copy conversion.

qPCR-scaled genus abundance was defined as genus relative abundance × `gene_copies_per_ul`. These values are estimated 16S target equivalents, not genus-specific qPCR assays or cell counts: 16S copy number varies within and between genomes and was not corrected here [Louca2018; Johnson2019], and relative enrichment does not imply absolute outgrowth [Props2017]. ⟦PLACEHOLDER: amplification efficiency, R², LOD/LOQ and qPCR technical variance were not reported in the provider Materials and Methods; obtain if needed for ALDEx3 `s.var`.⟧

**Outcomes**
Primary. 0-to-48-h changes and trajectories for acetate, propionate and butyrate.
zSecondary/supporting. Donor-aware 48-h endpoint concentrations for all five analytes; succinate and 5-aminovalerate trajectories; alpha diversity; Bray–Curtis community structure; total bacterial load.
Exploratory. Discovery-wide differential abundance; focal taxon–SCFA association; median-defined response strata.

**Statistical analysis**

Donor was the independent biological unit, with A/B aliquots nested within donor and culture wells nested within aliquot × carbohydrate linking paired 0- and 48-h observations [Lazic2010; Lazic2018]. Analyses requiring one observation per donor × condition averaged nested replicates; trajectory models retained the nested random structure below.

SCFA trajectories: Untransformed acetate, propionate and butyrate concentrations (mM) were fit with linear mixed-effects models of the form `concentration ~ group × carbohydrate × time + (1|donor) + (1|aliquot) + (1|well)` [Bates2015; Kuznetsova2017; Pinheiro2000; Zuur2009]. Type III tests used Satterthwaite degrees of freedom. Carbohydrate-specific 0-to-48-h changes and Tukey-adjusted pairwise differences in change (RDC vs no added carbohydrate, SDC vs no added carbohydrate, SDC vs RDC) were estimated with emmeans, marginalizing equally over obesity group [Lenth2017]. When residual non-normality or heteroscedasticity was evident, donor-bootstrap intervals accompanied model-based intervals; variance-component summaries also used bootstrap intervals because fixed-effect robustness under misspecification does not extend to random effects [Schielzeth2020].

SCFA concentrations at 48 h: After averaging nested replicates, endpoints were modeled as `concentration_48H ~ group × carbohydrate + (1|donor)`. Type III tests were BH-adjusted across 15 analyte × term tests [Benjamini1995]. These omnibus tests do not identify pairwise carbohydrate contrasts; plot-level Wilcoxon annotations were descriptive only.

Obesity-group trajectory contrasts: For each primary analyte under RDC and SDC, the estimand was `(obesity 48 h − obesity 0 h) − (healthy-weight 48 h − healthy-weight 0 h)`, with model-based Wald 95% and donor-bootstrap intervals. Equivalence testing was not performed because externally justified margins were unavailable.

Donor heterogeneity: Donor-level 0-to-48-h changes were modeled as `delta ~ group × carbohydrate + (1|donor)`. Repeatability across carbohydrate conditions was summarized by the ICC with donor-bootstrap intervals [Koo2016; Nakagawa2010]; RDC–SDC concordance by Spearman correlation with donor-bootstrap intervals; and obesity-group dispersion by SD ratios with bootstrap intervals and Brown–Forsythe tests BH-adjusted across six primary analyte × carbohydrate comparisons.

Diversity and community structure: Observed richness, Shannon and Simpson indices were averaged to donor × carbohydrate × timepoint and modeled as `metric ~ group × carbohydrate × time + (1|donor)`, with BH-adjusted obesity-group contrasts within carbohydrate × timepoint cells. Bray–Curtis dissimilarities on donor-aggregated relative abundances were tested by term-level PERMANOVA (9,999 permutations) [Anderson2001; McArdle2001; Anderson2017; Oksanen2001], with multivariate dispersion assessed by `betadisper` [Anderson2006; Anderson2013; Warton2012].

Differential abundance: Exploratory genus-level screens used FDR q < 0.1 and three complementary methods that were not combined into a single significance rule [Nearing2022; Nixon2025]. ANCOM-BC2 modeled relative abundances with a donor random intercept; carbohydrate and timepoint were encoded as a single `carb_time` factor, so coefficients are cell contrasts versus no added carbohydrate at 0 h [Lin2024]. MaAsLin 3 modeled log-transformed qPCR-scaled abundances (`relative abundance × gene_copies_per_ul`) as `~ carbohydrate × timepoint + obesity_group + (1|donor)` without additional normalization [Nickols2026]. ALDEx3 provided a scale-uncertainty sensitivity analysis with `s.mu = log2(gene_copies_per_ul)` and provisional `s.var ∈ {0.05, 0.25, 1.0}` (primary displays, `s.var = 0.25`) under the same fixed and random effects [McGovern2026].

**Integrated taxon–SCFA analysis.** The claim-facing association was the donor-level Spearman correlation between SDC-attributable net change in log10 qPCR-scaled *Fusicatenibacter* and SDC-attributable Δbutyrate after exact well matching and donor aggregation, with donor-permutation p-values and donor-bootstrap confidence intervals. Sample-level 48-h correlations and pooled well-level tests were retained only as descriptive diagnostics because they retain multiple records per donor.

**Response strata.** For exploratory visualization only, donors were assigned to within-analyte-and-carbohydrate median-split higher-response and lower-response strata after continuous analyses. Dichotomizing a continuous response discards information and yields unstable labels [Altman2006; Senn2018; Atkinson2015; Barnett2005]; strata were not interpreted as clinical or predictive phenotypes.

**Software.** Mixed models, diversity analyses and ANCOM-BC2 were run in R 4.3.3 (Bioconductor 3.18) with lme4 1.1.37, lmerTest 3.1.3, emmeans 1.11.2, vegan 2.7.1, phyloseq 1.46.0 and ANCOMBC 2.4.0 [Bates2015; Kuznetsova2017; Lenth2017; Oksanen2001; McMurdie2013; Lin2024]. MaAsLin 3 1.4.0 and ALDEx3 1.2.0 were fit in a separate locked R 4.6.0 / Bioconductor 3.23 project library [Nickols2026; McGovern2026]; focal taxon–SCFA and matching steps used Python. Versioned outputs and recorded `sessionInfo()` files are archived under `integrated/results/additional_analyses_2026-07-19/`; analysis scripts are in `integrated/integrated-scfa-microbiome.Rmd` and `integrated/additional_analyses/`. Code and result archives will be deposited upon publication (accession pending).

---

## 3. Results

### 3.1 Study sample

Sixteen independent donors (8 healthy-weight, 8 with obesity) contributed stool aliquots to the nested fermentation design (Figure 1; Table 1). After QC and matching, 289 SCFA–16S–qPCR records from these 16 donors entered claim-facing analyses (129 complete 0-to-48-h culture-well trajectories). Four of 20 planned donors were absent from the SCFA export ⟦PLACEHOLDER: reason⟧, and four design cells lacked a 48-h SCFA observation (including two obesity-group no-added-carbohydrate 48-h samples) ⟦PLACEHOLDER: reason and missingness assumption⟧. Donor × carbohydrate summaries and culture-well records are analytic units, not independent participants.

### 3.2 Endpoint SCFA concentrations differed by carbohydrate condition

In donor-averaged 48-h models, the omnibus carbohydrate term was BH-significant for propionate (F(2,28)=25.16, q=1.67×10⁻⁶), acetate (F(2,28)=16.11, q=6.63×10⁻⁵), butyrate (F(2,28)=8.58, q=0.0037) and succinate (F(2,28)=5.46, q=0.030), but not 5-aminovalerate (F(2,28)=1.41, q=0.78) (Figure 2; Table 2). No obesity-group main effect or group × carbohydrate interaction survived BH adjustment. These endpoint tests indicate that final concentrations differed among carbohydrate conditions and do not identify pairwise contrasts or 0-to-48-h accumulation differences.

### 3.3 SCFA trajectories and pairwise differences in change

Acetate, propionate and butyrate increased from 0 to 48 h under all three carbohydrate conditions, including no added carbohydrate (Figure 3A; Table 3). Mean within-condition changes were large for acetate (28.1–34.0 mM), moderate for butyrate (5.3–7.9 mM) and smaller for propionate (3.8–5.6 mM).

Adjusted pairwise differences in change were analyte-specific (Figure 3B; Table 3). Propionate accumulation was lower under RDC than under no added carbohydrate (−1.86 mM, 95% CI −3.57 to −0.15, adjusted p=0.025). The remaining eight contrasts were not significant at α=0.05, including SDC versus no added carbohydrate for butyrate (+2.62 mM, 95% CI −0.001 to 5.23, p=0.050) and acetate (+5.89 mM, 95% CI −0.23 to 12.01, p=0.066).

SDC did not differ significantly from RDC in 0-to-48-h change for acetate (+2.89 mM, 95% CI −2.06 to 7.84, p=0.54), propionate (+1.06 mM, 95% CI −0.32 to 2.44, p=0.23) or the primary butyrate contrast (+1.03 mM, 95% CI −1.08 to 3.15, p=0.73). Mean butyrate change was numerically higher under SDC (7.94 mM) than RDC (6.91 mM), but the adjusted difference was not significant.

### 3.4 Obesity status was not associated with SCFA change contrasts

None of the six primary obesity-group difference-in-change contrasts (acetate, propionate and butyrate under RDC and SDC) differed significantly from zero (Figure 4A; Table 4). Point estimates ranged from −1.19 to +3.43 mM, with the widest intervals for acetate. Donor-bootstrap intervals supported the same qualitative conclusion. Equivalence was not tested. Obesity-versus-healthy-weight dispersion ratios for donor-level changes were near one (Brown–Forsythe BH q ≥ 0.995).

### 3.5 Donor-level SCFA responses were repeatable across carbohydrates

Donor identity accounted for a substantial share of variation in 0-to-48-h SCFA change (Figure 4B; Table 5). ICCs across the three carbohydrate conditions were 0.84 (bootstrap 95% CI 0.49–0.95) for acetate, 0.84 (0.66–0.90) for propionate and 0.63 (0.41–0.76) for butyrate. Rank ordering of donors was highly concordant between RDC and SDC (Spearman ρ=0.96, 0.85 and 0.83 for acetate, propionate and butyrate, respectively). These patterns indicate generalized fermentative responsiveness rather than formulation-specific ranking.

### 3.6 Community structure tracked incubation time and obesity group more than carbohydrate condition

**Alpha diversity.** Donor-averaged Shannon, Observed and Simpson indices showed carbohydrate, time and carbohydrate × time associations in mixed models (Table 6). At 48 h under no added carbohydrate, obesity-group communities had lower Shannon (−0.50, 95% CI −0.84 to −0.17, p=0.004) and Observed richness (−42.8, 95% CI −69.8 to −15.8, p=0.002) than healthy-weight communities; Simpson did not differ (p=0.23). Corresponding group contrasts under RDC and SDC at 48 h were not significant.

**Beta diversity.** In donor-aggregated Bray–Curtis PERMANOVA across both timepoints, incubation time explained the largest share of variation (R²=0.239, p=0.0001), with a smaller obesity-group term (R²=0.026, p=0.0037); carbohydrate and group × carbohydrate were not significant (Figure 5; Table 7). At 48 h alone, group remained associated with community structure (R²=0.056, p=0.0027), whereas carbohydrate (R²=0.044, p=0.41) and the interaction (R²=0.010, p=1.0) were not. Multivariate dispersion at 48 h did not differ by carbohydrate (p=0.88) or group (p=0.094).

**Total bacterial load.** Donor-level Δlog₁₀ 16S gene copies increased more under RDC than under no added carbohydrate (+0.71, 95% CI 0.16 to 1.25, n=14, p=0.025). SDC versus no added carbohydrate (+0.37, 95% CI −0.08 to 0.83, p=0.13) and SDC versus RDC (−0.21, 95% CI −0.73 to 0.31, p=0.44) were not significant (Figure 5C; Table 8). Marked community expansion during incubation provides necessary context for relative-abundance contrasts.

### 3.7 Carbohydrate and time reshaped genus profiles on relative and qPCR-scaled scales

Exploratory genus screens (FDR q < 0.1) identified multiple taxa whose abundance changed with incubation and carbohydrate exposure; frameworks were treated as complementary estimands rather than a voting rule (Figure 6; Table 9).

**qPCR-scaled abundance (MaAsLin 3).** Five carbohydrate × 48-h interaction terms met the exploratory threshold, all positive: *Prevotella* under SDC (coefficient 18.09, joint q=1.7×10⁻⁵) and RDC (15.85, q=5.7×10⁻⁴); *Lactobacillus* under SDC (9.03, q=0.047); *Bifidobacterium* under SDC (3.28, q=0.085); and an unresolved Prevotellaceae feature under SDC (11.79, q=0.087). *Megasphaera* showed large main carbohydrate effects under both RDC and SDC (coefficients ≈20.2–20.4, joint q≈2×10⁻⁶). Seventeen genera also showed strong timepoint (48 h) main effects, including *Carnobacterium*, *Negativicoccus* and *Odoribacter*.

**Relative abundance (ANCOM-BC2).** Cell contrasts versus no added carbohydrate at 0 h identified overlapping 48-h enrichments under SDC and RDC for *Bacteroides*, *Parabacteroides*, *Clostridium*, *Escherichia*, *Finegoldia*, *Eggerthella* and *Phascolarctobacterium* (q < 0.05), with concurrent relative decreases in several Lachnospiraceae/Clostridia lineages including *Eubacterium*, *Anaerostipes*, *Dorea*, *Ruminococcus*, *Blautia* and *Marvinbryantia*. Obesity-group coefficients were not significant at q < 0.1. In a focused 48-h RDC-versus-SDC contrast, a single genus met q < 0.1 and was higher under SDC (*Fusicatenibacter*; Figure 6D) ⟦PLACEHOLDER: export head-to-head LFC, 95% CI and q⟧.

**Scale-uncertainty sensitivity (ALDEx3).** At provisional `s.var=0.25`, significant associations were confined to the 48-h timepoint main effect (e.g. *Bacteroides*, *Parabacteroides*, *Escherichia*, *Flavonifractor*, *Alistipes*, *Clostridium*, *Akkermansia*); no carbohydrate × 48-h interaction reached adj. p < 0.1. Directional SDC × 48-h estimates for *Fusicatenibacter* (+4.31) and *Bifidobacterium* (+3.26) were positive but nonsignificant after multiplicity adjustment.

**Cross-method concordance.** Among focus effects with all three coefficients available, 85/182 (47%) were three-way direction-concordant; for SDC × 48-h terms, 16/45 were three-way concordant, whereas MaAsLin 3 and ALDEx3 agreed in direction for 37/45. Several taxa (*Bacteroides*, *Parabacteroides*, *Clostridium*, *Escherichia*, *Eggerthella*, *Bifidobacterium*) were directionally positive on the qPCR-scaled interaction scale and also enriched in ANCOM-BC2 48-h cells, whereas some Lachnospiraceae showed opposing relative versus scaled signs—consistent with compositional change during community expansion rather than a single biological contradiction.

### 3.8 Donor-level *Fusicatenibacter*–butyrate association was positive but not significant

As a prespecified exploratory focal test, SDC-attributable net change in log₁₀ qPCR-scaled *Fusicatenibacter* was positively associated with SDC-attributable Δbutyrate at the donor level (Spearman ρ=0.48, permutation p=0.088, bootstrap 95% CI −0.24 to 0.89, n=14; Figure 7; Table 10). The interval included zero. Secondary specificity checks (RDC butyrate, relative *Fusicatenibacter*, total load, acetate, propionate and total SCFA) were not significant after BH adjustment within that family. Leave-one-donor-out Spearman estimates ranged approximately 0.35–0.63. Sample-level 48-h correlations retaining multiple wells per donor are reported only as descriptive diagnostics and are not used for inference.

### 3.9 Median-defined response strata were unstable exploratory labels

Within-analyte, within-carbohydrate median splits of donor-level SCFA change yielded higher- and lower-response labels that agreed across a donor’s paired A/B aliquots for only 7–12 of 16 donors (Supplementary Table S1). Community differences between butyrate or propionate strata were not significant by PERMANOVA (p=0.11 and 0.19, respectively). These data-derived strata are not interpreted as clinical phenotypes.

---

## 4. Discussion

### 4.1 Principal findings

Five findings frame this study. First, carbohydrate condition was associated with endpoint acetate, propionate, butyrate and succinate concentrations at 48 h in donor-aware models, with no BH-significant obesity-group or interaction term. Second, SCFAs accumulated during culture in all conditions, and carbohydrate condition modified trajectories in an analyte-specific manner, including lower propionate change under RDC than under no added carbohydrate; SDC did not differ significantly from RDC for any primary SCFA change, including butyrate. Third, obesity-group SCFA change contrasts were not significant, although equivalence was not established. Fourth, donor response ordering was strongly repeatable across the two active carbohydrates. Fifth, community structure tracked incubation time and, more modestly, obesity group, whereas carbohydrate condition was not associated with Bray–Curtis composition in claim-facing PERMANOVA.

Exploratory differential abundance implicated multiple genera rather than a single taxon: qPCR-scaled models highlighted *Prevotella*, *Lactobacillus*, *Bifidobacterium* and *Megasphaera* under carbohydrate exposure, while relative-abundance models showed concurrent 48-h enrichment of *Bacteroides*/*Parabacteroides* lineages and relative depletion of several Lachnospiraceae. A focused 48-h RDC-versus-SDC contrast singled out *Fusicatenibacter*, and the donor-level SDC-attributable *Fusicatenibacter*–butyrate association was positive with an interval including zero.

### 4.2 Relationship to the originating animal work

The rodent study that motivated this work incorporated carbohydrate formulations into complete diets, so host digestion determined colonic exposure and host metabolic outcomes were measurable [PlazaDiaz2022]. The present study bypassed digestion and exposed fecal communities directly. This is a difference in the independent variable, not only in the model system, and it has a specific consequence: the property that distinguishes a slow- from a rapid-digesting formulation — the rate at which host enzymes hydrolyze it — has no expression in a vessel containing no host enzymes. Isomaltulose, a major component of the slower-digesting blend, is almost completely digested and absorbed in the human small intestine [Holub2010]; static in vitro digestion cannot reproduce the kinetics that would matter here [Brodkorb2019]. A null SDC-versus-RDC contrast in this assay is therefore the expected result of the design rather than a failure to reproduce the animal finding.

Two further features of the anchor study bound the comparison. Its butyrate contrast was itself not significant, so the present butyrate null does not contradict it — a point that matters for how the prespecified butyrate hypothesis is justified.

More importantly, its slower-digesting arm varied several things at once. The cellulose fraction was exchanged for resistant maltodextrin plus inulin/FOS, and the sugar base was changed from sucrose plus cornstarch to isomaltulose plus sucromalt, with higher total sugar in the slower-digesting diet (39.10 vs 36.85 g/100 g) and an approximately three-fold difference in glycemic load. Digestion rate, fermentable-substrate identity and sugar composition therefore covary, and no contrast in that design isolates digestion rate. ⟦Note for authors: the anchor's Table 1 also contains an internal arithmetic inconsistency — the declared fiber components sum to 17.00 rather than the stated 16.00 g/100 g. Verify against the published table before citing any composition figure.⟧

The animal work is best read as biological rationale and a hypothesis generator, not as a replication target [Walter2020].

### 4.3 Interpreting carbohydrate effects

The endpoint and change estimands answer different questions and should not be merged. The 48-h omnibus tests establish that the final metabolic state varied somewhere among carbohydrate conditions for four analytes; they do not identify which condition was higher, nor whether any difference arose from greater 0-to-48-h accumulation.

All conditions accumulated SCFAs, including the no-added-carbohydrate control. This is expected rather than anomalous. Peptone-based fermentation media and the fecal inoculum are themselves substantial nutrient sources [PerezBurillo2021], and "no added carbohydrate" does not mean carbohydrate-free: batch-fermentation basal media of this class routinely contain starch, peptone, yeast extract and mucin at gram-per-litre concentrations, all of which are fermentable [VandenAbbeele2022]. Non-carbohydrate substrates can also be fermented to SCFAs — in a piglet colonic model, glutamine and glutamate each raised acetate by roughly 20 mM and butyrate by roughly 10 mM above a no-substrate blank, whereas branched-chain amino acids *lowered* acetate [VandenAbbeele2022]. A rich background can additionally compress the measurable effect of an added substrate: reducing basal medium nutrient density increases the apparent impact of added fiber [Poppe2023].

The one significant pairwise change contrast — lower propionate under RDC than under no added carbohydrate — argues directly against a "more substrate equals more SCFA" narrative. Several mechanisms could produce an analyte-specific or negative contrast, and we present them as hypotheses rather than conclusions. Because medium pH was set at the start of incubation but not maintained, acid accumulation over 48 h would be expected to lower pH differentially across conditions, and falling pH suppresses the succinate-pathway propionate producers among *Bacteroides* while favoring butyrogenesis [Walker2005; Duncan2009; LaBouyer2022]. A readily fermented substrate can also suspend cross-feeding until it is exhausted, as glucose does for lactate utilization [Duncan2004], and end-product inhibition can depress multiple SCFAs together [Wang2020]. Propionate pathways are, empirically, often the least responsive to added fiber [VanWehle2024]. We did not measure pH, substrate disappearance or flux, so none of these mechanisms is tested here; distinguishing them is a specific aim for follow-up.

Finally, a single 48-h endpoint may under-resolve the very property of interest. Substrates with similar endpoint SCFA can follow entirely different fermentation time courses within a 48-h window [Kaur2011; Rose2009; Gu2018]. A null endpoint contrast is therefore compatible with a real difference in fermentation *rate* that this design cannot see.

### 4.4 Obesity-group results

None of the six obesity-group difference-in-change contrasts was significant. This is a **non-detection, not a demonstration of equivalence**: with 8 donors per group and no externally justified equivalence margin, the study is not able to distinguish "no meaningful difference" from "a difference this design could not resolve." The acetate intervals in particular remain wide. We therefore avoid describing fermentation capacity as preserved, intact, comparable or equivalent, and note that a claim in that direction would require a prespecified equivalence test with margins justified independently of the observed contrasts.

The surrounding literature does not point in a single direction either. Adult observational data lean toward *higher* fecal SCFA in obesity, with substantial heterogeneity [Kim2019]; pediatric studies split in all three directions; and ex vivo work has found lean-versus-obese differences that reverse sign depending on the substrate tested [Aguirre2014]. Our nonsignificant contrasts are thus unsurprising without being informative.

A caution from the closest methodological predecessor is worth carrying forward: in adolescents with obesity, fecal SCFA concentration and ex vivo SCFA production capacity were *negatively* associated, and neither tracked obesity markers in the expected direction [Holmes2020]. Fecal concentration and fermentation capacity are therefore not interchangeable readouts, which is a further reason not to attach a health valence to any observed SCFA difference in either direction.

It is worth stating plainly that the community and metabolite results give different group signals: obesity group was associated with Bray–Curtis structure, while SCFA change contrasts were not significant. Composition and measured function need not track one another [Reichardt2018; Oliver2021; Yao2024]. Functional redundancy is one candidate explanation, but it should not be asserted here — the statistical signature of redundancy can arise from averaging alone [Ho2025], and invoking it would amount to the equivalence claim we have declined to make.

### 4.5 Interpersonal response heterogeneity

The most robust positive finding in this study is the magnitude and repeatability of between-donor variation. Donor ICCs of 0.63–0.84 and RDC–SDC rank correlations of 0.83–0.96 indicate that donors occupied consistent positions in the high-to-low response ordering irrespective of which active formulation they received. This is best described as **generalized fermentative responsiveness** rather than substrate-specific response.

Several alternatives are not excluded by these data and should be conceded alongside the estimates rather than deferred. Shared donor-level nuisance scaling — inoculum density, baseline SCFA pool — would produce cross-condition correlation without any shared biology of substrate use. RDC and SDC may simply be similar enough as microbial substrates that a high correlation is uninformative about response specificity. And because each donor contributed a single fermentation occasion, donor variance and donor-by-run variance are not separable here.

We note explicitly that this result does **not** validate a responder phenotype, a classifier, or a dietary-response predictor. Median-split labels applied to a donor's own paired aliquots agreed for only 7–12 of 16 donors, which is consistent with a general property of dichotomized continuous responses: competing classification rules applied to identical data classify a substantial minority of individuals inconsistently [Hecksteden2018; Atkinson2015; Senn2018]. The appropriate framing for future work is a continuous phenotype with repeated measurement, not a binary label.

A tension in the literature deserves acknowledgment. In one in vitro study, "despite marked inter-individual differences in microbiota composition, SCFA production was surprisingly reproducible for different carbohydrates," implying that functional redundancy constrains between-person variation [Reichardt2018] — though that observation rests on only three donors, which limits how much weight it can bear against the present 16-donor estimates. The reconciliation most consistent with available evidence is analyte-specific: redundancy appears to dominate for total SCFA and acetate, while propionate and butyrate vary more substantially between individuals [Fan2026]. Our own ICCs are consistent with that pattern, being lowest for butyrate.

### 4.6 *Fusicatenibacter* as a hypothesis, not a mechanism

An exploratory relative-abundance analysis distinguished SDC from RDC at 48 h by a single genus. Several considerations constrain what may be concluded.

The direction of the finding is consistent and is worth stating plainly: on the load-aware scale, *Fusicatenibacter* increased during incubation, and increased more under SDC than under either the paired control or RDC. Because total bacterial load itself rose 5- to 24-fold across conditions, this is a case where relative abundance is a poor guide to what happened — the scaling to total load was necessary rather than decorative, and it is the reason the study measured load at all [Props2017].

What the data do not support is a claim of statistical significance. The multivariable qPCR-scaled coefficients were substantial but did not survive multiplicity adjustment (MaAsLin 3 q=0.125; ALDEx3 q=0.581), which with 14–15 donors and wide between-donor spread reflects imprecision rather than absence of effect. The donor-level association between SDC-attributable *Fusicatenibacter* change and butyrate change was positive but its interval included zero. The relative- and load-aware analyses agree in direction once the ANCOM-BC2 cell parameterization is read against the paired control (§3.8), so the frameworks are consistent — but consistency of direction across screens is not a significance test, and these remain exploratory genus-wide analyses.

A useful specificity check argues against the simplest artifactual explanation: the association between SDC-attributable butyrate change and total bacterial load change was slightly negative (ρ=−0.121), so the *Fusicatenibacter*–butyrate relationship is not merely a reflection of communities that grew more overall producing more of everything. Changes in taxon abundance and metabolite concentration were measured concurrently, so direction of causation is not established. Genus-level 16S data do not establish strain-level function, and butyrate synthesis is polyphyletic and not inferable from 16S taxonomy [Vital2014]. Substrate disappearance and metabolite flux were not measured.

One point requires more than a hedge. *Fusicatenibacter* is **not a described butyrate producer**: the type-strain description of *F. saccharivorans* reports lactic, formic, acetic and succinic acid as fermentation end products from glucose, and butyrate is not among them [Takada2013]. We are not aware of isolate, genomic or flux evidence establishing butyrogenic capacity in this genus. It would therefore be incorrect to interpret this signal as identifying the source of butyrate in these cultures. A more defensible hypothesis is **cross-feeding donation** — that *Fusicatenibacter* fermentation products, particularly lactate and acetate, are substrates for butyrogenic partners such as *Anaerostipes* or *Eubacterium hallii* [Duncan2004; Belenguer2006] — but this was not tested here and is offered only as a candidate mechanism.

Context makes the observed pattern unsurprising. *Fusicatenibacter* is carbohydrate-responsive, increasing with resistant maltodextrin supplementation [Mai2022], and its abundance has shifted without concurrent SCFA change in independent studies. A relative signal lacking absolute-scale corroboration is the previously observed norm for this genus rather than an anomaly requiring mechanistic explanation. Its human metabolic associations are inconsistent in direction; in the closest demographic comparison, adolescents with obesity, it was *higher* in those with NAFLD [Testerman2022].

Appropriate follow-up would include isolate growth on the exact RDC and SDC formulations, substrate disappearance measurement, metabolite profiling of pure and co-culture, and ideally ¹³C tracing of the proposed cross-feeding route [Belenguer2006; Ze2012].

### 4.7 Strengths

Each donor's community was exposed to all three conditions in a paired design, so carbohydrate comparisons are within-donor. The analysis respected the donor → aliquot → culture-well hierarchy and treated donor as the biological unit, avoiding the pseudoreplication that inflates precision in nested fermentation designs [Lazic2018]. Targeted metabolite quantification, 16S profiling and total bacterial qPCR were integrated, and differential abundance was examined on both relative and qPCR-scaled scales as distinct estimands. Contrasts are reported as adjusted estimates with uncertainty rather than as significance verdicts alone, and donor repeatability was quantified rather than asserted. The donors were adolescents, an age group for which ex vivo carbohydrate fermentation data are notably scarce.

### 4.8 Limitations

The ex vivo system lacks host digestion, absorption, transit, mucosal surfaces, immune signaling and dietary context; results speak to microbial response potential, not to in vivo delivery or clinical benefit. Expected host digestibility is a property of the formulation and was not measured in culture. SCFA concentrations reflect net accumulation rather than production flux [Sakata2019]. Only 0- and 48-h timepoints were analyzed, and 48 h entails substantial culture adaptation as well as substrate use; a single endpoint may under-resolve rate differences [Kaur2011]. Medium pH was set at the outset but not maintained, so pH drift is an uncontrolled and unmeasured covariate that is partially entangled with carbohydrate condition in its effect on butyrate and propionate stoichiometry [Duncan2009]. The cohort comprised 16 donors, 8 per group, with unresolved exclusions from the planned 20. Cryopreserved inocula are an accepted but not equivalent substitute for fresh material, with known effects on quantitative SCFA output [Aguirre2015]. SCFA concentration dilution mapping, CV and LOQ require facility confirmation. 5-aminovalerate may rest on a different measurement basis from the four calibrated analytes, which is an alternative explanation for its null carbohydrate term that has not been excluded. Total bacterial qPCR technical variance was unavailable, so the scale-uncertainty sensitivity analysis used an assumed rather than measured variance grid. Residuals were non-normal for all three primary analytes and heteroscedastic for propionate and butyrate. Genus-level 16S carries compositional and taxonomic-resolution limits, and qPCR-scaled estimates inherit uncertainty from both relative abundance and total-load measurement without 16S copy-number correction [Louca2018]. Agreement among differential-abundance frameworks was limited, partly because they estimate different quantities. The differential-abundance models accounted for donor but not the A/B aliquot split, because the dual random-effect specification failed to fit; the SCFA models did include it, so the two analysis families do not share identical nesting. Taxon and metabolite changes were concurrent, limiting causal and predictive interpretation. Exploratory analyses used a q<0.1 threshold across multiple analysis families. Median dichotomization produced labels sensitive to aliquot variation. No in vivo or clinical metabolic outcome was measured.

### 4.9 Conclusion

At 48 hours, carbohydrate condition was associated with endpoint acetate, propionate, butyrate and succinate concentrations in donor-aware models. Adjusted 0-to-48-hour changes under SDC did not differ significantly from RDC for acetate, propionate or butyrate, and obesity status did not significantly modify SCFA changes in this cohort, although equivalence was not established. Ex vivo accumulation also varied by incubation time and, substantially and repeatably, by donor identity. Community and taxon analyses generated testable hypotheses; mechanism and in vivo relevance require independent validation.

---

## Figure legends

### Figure 1. Study design and analytic sample flow

Schematic of the paired ex vivo fermentation design and analytic accounting. Sixteen adolescent donors (8 healthy-weight, 8 with obesity) each contributed one stool sample, divided into paired A and B aliquots. Each aliquot was incubated under three conditions — no added carbohydrate, RDC and SDC — in replicate culture wells (R1/R2, S1/S2), sampled at 0 and 48 h. The diagram distinguishes biological participants (donors) from nested technical and culture units (aliquots, wells); aliquots and wells are not independent participants. Counts shown: 16 donors, 32 aliquots, 160 culture wells in the nested SCFA workflow, 292 microbiome profiles matched to qPCR records. Exclusions and missingness are annotated, including four of 20 planned donors absent from the SCFA export and two missing obesity-group no-added-carbohydrate 48-h observations. ⟦PLACEHOLDER: finalize exclusion reasons.⟧

### Figure 2. Donor-aware SCFA concentrations at 48 hours

Donor-level distributions of acetate, propionate, butyrate, succinate and 5-aminovalerate concentrations at 48 h by carbohydrate condition, after averaging A/B aliquots and culture wells to one value per donor × carbohydrate × analyte. Points are donors; boxes show median and interquartile range. Panels are annotated with the omnibus Type III carbohydrate F statistic and BH-adjusted q-value from the model `concentration_48H ~ group × carbohydrate + (1|donor)`. **This figure does not resolve pairwise direction**: the omnibus test indicates that endpoint concentrations differed somewhere among conditions and does not identify which pair differed. Descriptive pairwise Wilcoxon annotations are deliberately omitted.

### Figure 3. SCFA trajectories and adjusted carbohydrate contrasts

**(A)** Estimated marginal trajectories of acetate, propionate and butyrate from 0 to 48 h by carbohydrate condition, with 95% confidence bands, from the nested mixed model `concentration ~ group × carbohydrate × time + (1|donor) + (1|aliquot) + (1|well)`. Accumulation in the no-added-carbohydrate condition is plotted on the same axes to make background fermentation visible. **(B)** Forest panel of the nine Tukey-adjusted pairwise differences in 0-to-48-h change (three contrasts × three analytes), with 95% confidence intervals and a reference line at zero. Note the single significant contrast — propionate lower under RDC than under no added carbohydrate — and that the SDC-versus-control butyrate interval crosses zero (p=0.0502).

### Figure 4. Obesity-group contrasts and donor repeatability

**(A)** Forest plot of the six obesity-minus-healthy-weight differences in 0-to-48-h change (acetate, propionate, butyrate × RDC, SDC), showing model-based and donor-bootstrap 95% intervals with a reference line at zero. None differs significantly from zero; equivalence was not evaluated, and the intervals should not be read as establishing it. **(B)** Donor-level scatterplots of RDC versus SDC 0-to-48-h change for each primary analyte, one point per donor, annotated with the donor ICC and RDC–SDC Spearman ρ (both with bootstrap 95% intervals). The tight positive association indicates repeatable generalized response ordering across formulations, not preferential response to either formulation.

### Figure 5. Community structure

Principal coordinates analysis of Bray–Curtis dissimilarities on donor-aggregated data. **(A)** Samples colored by incubation time (0 vs 48 h). **(B)** Samples colored by obesity group. Term-specific PERMANOVA R² and p-values are annotated on each panel (time R²=0.239, p=0.0001; group R²=0.026, p=0.0037 across all timepoints; group R²=0.055, p=0.0027 at 48 h). **(C)** Optional 48-h panel colored by carbohydrate condition — **the caption must state that the carbohydrate term was not significant at 48 h (R²=0.0435, p=0.409)**, and this panel should not be read as showing a carbohydrate effect. Multivariate dispersion tests at 48 h: carbohydrate p=0.876, group p=0.094.

### Figure 6. Exploratory *Fusicatenibacter* triangulation

**Labeled exploratory.** **(A)** Focused ANCOM-BC2 estimate for relative *Fusicatenibacter* abundance, SDC versus RDC at 48 h, with 95% interval and the exploratory q<0.1 threshold marked ⟦PLACEHOLDER: exact coefficient/interval/q pending versioned export⟧. **(B)** Total bacterial 16S load change from 0 to 48 h by condition (mean Δlog₁₀ +0.708 control, +1.176 SDC, +1.385 RDC), plotted alongside donor-level net change in qPCR-scaled *Fusicatenibacter* (+1.096 SDC, +0.837 RDC). This panel is essential context: the community expanded 5- to 24-fold, so relative abundance alone cannot show what happened to this genus, and the larger SDC *Fusicatenibacter* increase occurred despite a *smaller* net total-load increase than RDC. **(C)** qPCR-scaled SDC×48-h coefficients from MaAsLin 3 (3.777, q=0.125) and ALDEx3 at `s.var=0.25` (4.308, q=0.581) with intervals — both positive and substantial, neither significant after multiplicity adjustment. If the ANCOM-BC2 relative estimates are also plotted, difference the 48-h cells against the paired 48-h control rather than showing raw coefficients against the 0-h reference cell — the raw form reads as negative and invites the false conclusion that SDC reduced the genus (§3.8). **(D)** Donor-level scatterplot of SDC-attributable net change in qPCR-scaled *Fusicatenibacter* against SDC-attributable net change in butyrate (N=14 exactly well-matched donors), annotated with Spearman ρ=0.477, permutation p=0.088 and bootstrap 95% CI −0.239 to 0.892.

The caption must state that *Fusicatenibacter* is **not a described butyrate producer** — the type-strain description reports lactic, formic, acetic and succinic acid as end products from glucose [Takada2013] — and that no substrate-use, isolate-growth or flux measurement was performed. Consider moving the whole figure to the Supplement if the focused ANCOM estimate cannot be exported reproducibly.

---

## Table placeholders

### Table 1. Participant characteristics and assay accounting

⟦PLACEHOLDER — populate from clinical and facility records.⟧ Donor-level demographics by obesity group (n=8 per group): age, sex, race/ethnicity, BMI and BMI percentile, obesity definition and reference standard, pubertal stage if collected, antibiotic exposure history, chronic medications, dietary restrictions. Assay accounting: donors, aliquots, culture wells and microbiome profiles by condition and timepoint; missingness and exclusions with reasons. Footnote must state that A/B aliquots and culture wells are nested technical units, not participants.

### Table 2. Primary SCFA estimates

Four labeled panels, kept separate because they answer different questions:

- **Panel A — 48-h endpoint omnibus tests.** All 15 donor-aware Type III analyte × term tests (carbohydrate, group, group×carbohydrate for five analytes) with F, df, raw p and BH q. Source: `integrated/results/additional_analyses_2026-07-19/scfa_48h_endpoint_type3_anova.csv`.
- **Panel B — within-condition 0-to-48-h changes** for each analyte and condition, with 95% intervals and analyzed donor count.
- **Panel C — adjusted pairwise carbohydrate differences in change**, nine rows, with 95% intervals and Tukey-adjusted p-values. ⟦Refresh acetate and propionate rows from the repaired live nested model before locking.⟧
- **Panel D — obesity-group difference-in-change contrasts**, six rows, with model-based and donor-bootstrap 95% intervals. Footnote: equivalence not evaluated.

Also report donor ICCs and RDC–SDC Spearman correlations with bootstrap intervals (or place in a separate Table 2E). Footnote that concentrations are reported in mM.

### Table 3. Community and focal-taxon inference

⟦PLACEHOLDER — assemble from integrated analysis exports.⟧ Alpha-diversity contrasts with the locked multiplicity family; term-specific PERMANOVA R² and p-values for all-timepoint and 48-h donor-aggregated models; betadisper results; focused ANCOM-BC2 estimate for *Fusicatenibacter* (SDC vs RDC, 48 h) with interval and q; MaAsLin 3 and ALDEx3 qPCR-scaled SDC×48-h estimates; donor-level *Fusicatenibacter*–butyrate association (ρ, permutation p, bootstrap CI, N). Footnote that frameworks estimate different quantities and that concordance was not used as a significance criterion.

---

## References

Full bibliographic entries — with PMIDs and DOIs, verification status, and claim-lock restrictions — are in **`manuscript/references.md`**, compiled from the 18 verified research reports in `research/`.

**Before finalizing, note three version traps** documented in that file: cite MaAsLin 3 as Nickols et al., *Nature Methods* 2026 (not the superseded bioRxiv preprint); cite ALDEx3 via McGovern & Silverman, *Microbiome* 2026 plus the software artifact (it has no standalone methods paper); and cite ANCOM-BC2 as the *Nature Methods* paper, not the Research Square preprint.

**Citations restricted by the claim lock** (see `research/INDEX.md` §4) include Li 2015, Bartsch 2025, Gerasimidis 2022, Reichardt 2018, Holmes 2020/2022, Gurry 2021 and Zeevi 2015. Each may be cited only for the narrow purpose noted there — none may be cited in support of SDC superiority, obesity equivalence, preserved fermentation capacity, responder phenotypes, or dietary-response prediction.

All 88 in-text citation keys in this draft have been cross-checked and resolve to entries in `references.md`. Two notes: `Reichardt2017` and `Reichardt2018` are the **same** paper (ISME J, PMID 29192904), and only `Reichardt2018` is used here; and `Haak2018` was added to the bibliography from the DFI-HMMF methods document rather than from the literature reports, so its identifier has not been independently re-verified.

---

## Outstanding items before this draft can advance

Carried forward from `manuscript/manuscript-outline.md` §"Pre-drafting checklist" and updated with what this draft surfaced:

1. Reconcile the four planned donors absent from the SCFA export and the two missing obesity no-carb 48-h observations.
2. Confirm participant characteristics, obesity definition, ethics approvals and exclusions.
3. Obtain the full fermentation protocol: MIPRO formulation and starting pH, RDC/SDC composition and dose, inoculum loading, atmosphere, temperature, vessel volume, agitation, termination.
4. **Concentration unit resolved: mM.** The earlier micromolar label was incorrect and is being corrected in the integrated analysis. Verify that all downstream tables, figures and exports have been regenerated, including the `margin_lower_uM` / `margin_upper_uM` column names in `scfa_metabolomics/equivalence_margins.csv` and the `unit_status` field of the equivalence-status export. Still outstanding from the facility: project-specific %CV, internal-standard recovery and LOD/LOQ; state the below-LOQ rule and confirm it was applied uniformly across analytes and arms.
4a. **Confirm whether 5-aminovalerate was calibrated** or reported as normalized relative abundance. The facility calibration curve covers only acetate, propionate, butyrate and succinate — and 5-aminovalerate is the single analyte with a null carbohydrate term, so a differing measurement basis is a live alternative explanation.
4b. **Confirm the GC-MS instrument.** The facility overview says Agilent 8890; the detailed PFBBr method says 7890A/5975C. The draft uses the latter.
4c. **Request qPCR technical variance from Zymo** so the ALDEx3 scale-uncertainty analysis uses a measured rather than assumed `s.var`.
4d. **State which taxonomy assignment** (strobelca, sintax or raxtax against SILVA 138.2) was carried into the claim-facing analysis, and whether the three agree at genus level for *Fusicatenibacter*.
5. Refresh the versioned SCFA pairwise export from the repaired live nested model; reconcile the remaining butyrate contrasts.
6. Repair and re-render the failed integrated ANCOM sections; export the focused 48-h RDC-versus-SDC ANCOM-BC2 coefficient, interval and q-value.
7. Lock the alpha-diversity multiplicity family.
8. Confirm betadisper was run at 48 h (reported here) and report it alongside every PERMANOVA term.
9. Reconsider the stated rationale for treating SDC-versus-RDC butyrate as a primary hypothesis: the rodent anchor's own butyrate contrast was null, so prespecification — not the anchor paper — must carry that choice.
10. Audit all draft text for any implication that *Fusicatenibacter* **produces butyrate** — the type-strain description reports lactic, formic, acetic and succinic acid as end products, not butyrate. This restriction is about metabolic function only; the abundance finding itself is directionally positive under SDC on the load-aware scale and should be reported as such.
10b. When reporting ANCOM-BC2, **state the `carb_time` cell parameterization explicitly** and difference the 48-h cells against the paired control before comparing with MaAsLin 3 / ALDEx3 interaction terms. Read as raw coefficients against the 0-h reference cell they look negative and invite the false conclusion that SDC reduced the genus; differenced against the 48-h control they agree with every other analysis. The `absolute_da/README.md` already documents this — carry that note into the manuscript Methods so a reader cannot repeat the error.
10a. Annotate the two `*_ERROR.txt` artifacts in `integrated/results/additional_analyses_2026-07-19/absolute_da/`, which mean different things. `maaslin3_ERROR.txt` is genuinely stale — a later run completed successfully and wrote the results. `ancombc2_ERROR.txt` is **not** stale: it records the interaction-parse failure that led to the `carb_time` combined-factor workaround actually used, and so is part of the methods provenance. Label each accordingly so neither is misread.
11. Position the closest human precedent accurately — Gillen et al. 2021, *Clinical Nutrition* (PMID 34130017, NCT03185884), Abbott-affiliated. It is **published**, not a preprint ("pre pub" in the local filename means *pre-pubertal*). Three details to state correctly: its SDC blend is essentially identical to the anchor study's (isomaltulose 26.0 vs 26.40, sucromalt 22.0 vs 22.10, resistant maltodextrin 10.0 vs 10.00, inulin-FOS 7.0 vs 7.00; Gillen adds isomaltooligosaccharide and reduces maltodextrin), but the **rapid-digesting comparators differ** (a 50:50 maltodextrin:sucrose drink versus a rodent high-fat diet); it measured substrate oxidation **plus blood glucose, insulin and DXA body composition**, so the accurate gap statement is that it has **no microbiome or SCFA endpoint**, not that it measured only oxidation; and inclusion was BMI >5th and <95th percentile, which admits overweight children, so describe the cohort as **non-obese**, not healthy-weight.
12. Complete author list, affiliations, funding and conflict-of-interest disclosures.
