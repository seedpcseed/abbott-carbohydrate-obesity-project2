# Deep-research citation prompts

**Project:** Abbott–Lurie ZC07, Project 2  
**Date:** July 19, 2026  
**Evidence sources:** `manuscript/manuscript-outline.md`, `integrated/integrated-scfa-microbiome.html`  
**Purpose:** Generate leading citations and prior claims for Introduction, Methods, and Discussion drafting—without expanding unsupported manuscript claims.

## How to use these prompts

Paste each prompt into a deep-research LLM (Gemini Deep Research, ChatGPT Deep Research, Perplexity, Consensus, or similar). Prefer **one prompt per research run**. After each run, keep:

1. **5–12 primary citations** (prefer reviews + key primary papers; peer-reviewed; PMID/DOI)
2. **Prior claims** stated as short propositions with supporting citations
3. **Conflicts / nuance** where the literature disagrees with a naive reading of our results
4. **What we must not overclaim**, given our claim lock

### Shared study facts to append when useful

Use these if a prompt asks for study-specific context:

> Human adolescent ex vivo fecal fermentation study (n=16 donors; 8 healthy-weight, 8 with obesity). Each donor’s fecal microbiota was incubated with no added carbohydrate, a rapid-digesting carbohydrate formulation (RDC), or a slow-digesting carbohydrate formulation (SDC). Primary endpoints were acetate, propionate, and butyrate net accumulation (0→48 h). 16S V3–V4 community profiling and total bacterial 16S qPCR scaling were secondary/exploratory. This is functional phenotyping of microbial communities, not an in vivo efficacy, digestion, obesity-treatment, or host-glycemic study. Originating rationale: Plaza-Díaz et al. 2022 rodent dietary carbohydrate-replacement study.

### Global language guardrails for all runs

Do **not** retrieve or invent citations that would support: SDC superiority over RDC for SCFAs; obesity equivalence / “preserved fermentation capacity”; absolute *Fusicatenibacter* expansion; *Fusicatenibacter* as a proven SDC utilizer or butyrate producer; median “responder phenotypes”; microbiome prediction of dietary response; or clinical benefit of RDC/SDC from this design.

---

## Prompt package (18)

### Introduction

#### Prompt 01 — Pediatric obesity and metabolic-response heterogeneity

**Manuscript target:** Introduction ¶1  
**Citation job:** Clinical framing without overselling microbiome causality.

```text
Research prompt for literature review and citation assembly.

TOPIC
Pediatric / adolescent obesity, with emphasis on interpersonal heterogeneity in metabolic and dietary responses—not a catalog of obesity-associated microbiome signatures.

NEEDED OUTPUT
1. Concise background claim set (8–12 short propositions) suitable for one Introduction paragraph.
2. Leading citations for:
   - Pediatric/adolescent obesity prevalence and metabolic comorbidity burden
   - Heterogeneity of weight-loss, glycemic, and diet-response outcomes in youth
   - Why between-person variation matters for nutrition and microbiome-targeted interventions
3. Explicitly distinguish host metabolic phenotype from microbial community phenotype.
4. Prefer recent reviews (2018–2026) plus seminal primary papers; include DOI/PMID.
5. Flag contested claims (e.g., microbiome as cause vs correlate of pediatric obesity).

CONSTRAINTS
Do not center the paragraph on “dysbiosis lists.” Do not imply our study treats obesity or measures host glycemic benefit.
Return: prioritized citation table (title, year, journal, DOI/PMID, why cite), plus 8–12 claim bullets each mapped to citations.
```

#### Prompt 02 — Microbial carbohydrate fermentation and SCFA net accumulation

**Manuscript target:** Introduction ¶2; Discussion 4.3  
**Citation job:** Correct biology of acetate/propionate/butyrate; net accumulation vs production/flux.

```text
Research prompt for literature review and citation assembly.

TOPIC
Gut microbial fermentation of nondigested carbohydrates to short-chain fatty acids (acetate, propionate, butyrate), emphasizing that measured culture or fecal concentrations often reflect net accumulation (production minus consumption/conversion), not metabolic flux or host systemic exposure.

NEEDED OUTPUT
1. Core physiology citations: major producing taxa and pathways; cross-feeding; butyrate as secondary fermentation product.
2. Citations warning against equating SCFA concentration with “production,” glycemic benefit, or colonic delivery.
3. Key reviews and primary papers on SCFA roles in gut physiology and host metabolism (balanced: potential benefits AND caveats).
4. Claim set distinguishing:
   - substrate fermentation potential
   - net SCFA accumulation in culture
   - fecal SCFA in vivo
   - circulating SCFA / host effects
5. Note typical time scales of batch fermentation assays.

CONSTRAINTS
Our manuscript must say “net accumulation / net output,” not “production,” unless carefully defined. Prefer human-relevant or methodologically rigorous sources. Include DOI/PMID.
Return: claim bullets + prioritized citations + 3–5 “language pitfalls to avoid.”
```

#### Prompt 03 — Carbohydrate structure, digestibility class, and microbial fermentation outcomes

**Manuscript target:** Introduction ¶3; Discussion 4.2–4.3  
**Citation job:** Why rapid- vs slow-digesting carbohydrate formulations can matter microbiologically—without claiming our head-to-head superiority.

```text
Research prompt for literature review and citation assembly.

TOPIC
How carbohydrate molecular structure, digestibility kinetics, and processing affect gut microbial fermentation endpoints (SCFA profiles, taxon shifts), contrasting:
A) in vivo dietary systems where host digestion filters colonic substrate availability
B) ex vivo / in vitro fecal fermentation systems where substrate is applied directly to microbiota

NEEDED OUTPUT
1. Leading citations on resistant starch, slowly digestible starch, and multi-component carbohydrate blends vs rapidly digestible starch regarding microbial SCFA responses.
2. Method papers / reviews contrasting batch culture vs SHIME / TIM-2 / continuous systems.
3. Explicit discussion of why “slow-digesting” as a host-digestibility label may not map 1:1 onto ex vivo microbial exposure.
4. Claim set explaining expected analyte-specific SCFA responses (acetate vs propionate vs butyrate) rather than uniform SCFA increases.
5. Identify authoritative reviews that caution against “more substrate → more SCFA” as a universal rule (background medium, pH, utilization, cross-feeding).

ANCHOR PAPER TO CONTEXTUALIZE (do not treat as our result)
Plaza-Díaz et al., Front Nutr 2022 (DOI 10.3389/fnut.2022.992682): rodent high-fat diet with rapid vs slower-digesting carbohydrate replacement; microbiome and fecal acetate/propionate changes.

CONSTRAINTS
Do not retrieve literature framed as proving SDC superiority for human adolescents. Separate host digestion claims from microbial fermentation claims.
Return: claim bullets, citation table, and a short “design translation caveats” paragraph suitable for Discussion 4.2.
```

#### Prompt 04 — Ex vivo / in vitro human fecal fermentation as a phenotype assay

**Manuscript target:** Introduction ¶4–5; Discussion 4.7–4.8  
**Citation job:** Justify the assay class and its inferential limits.

```text
Research prompt for literature review and citation assembly.

TOPIC
Use of human fecal microbiota in anaerobic ex vivo / in vitro batch fermentation to phenotype interpersonal variation in carbohydrate fermentation and SCFA accumulation.

NEEDED OUTPUT
1. Canonical methods/reviews for fecal-slurry batch fermentation (medium choices, inoculum loading, anaerobic conditions, sampling times).
2. Prior human studies—ideally including adolescents or children if available—showing that donors differ substantially in SCFA response to the same carbohydrate challenge.
3. Citations establishing what such assays can and cannot support:
   CAN: relative microbial response potential; donor heterogeneity; formulation comparisons under direct exposure
   CANNOT: host digestion/absorption; glycemic outcomes; mucosal/immune interactions; clinical efficacy
4. Discussion of culture-time adaptation (why 24–48 h endpoints reflect community adaptation as well as substrate use).
5. Preferred citations for “functional phenotyping” framing of microbiome studies.

CONSTRAINTS
Prefer human over animal method papers when available. Flag pediatric gaps. Include DOI/PMID.
Return: 10 prioritized methods/phenotype citations + claim bullets for translational-gap paragraph.
```

#### Prompt 05 — Obesity status and fecal fermentation capacity: prior evidence and null/ambiguous results

**Manuscript target:** Introduction ¶5; Discussion 4.4  
**Citation job:** Prior claims that obesity impairs fermentation—and contrary or null human evidence—to frame our nonsignificant group contrasts without claiming equivalence.

```text
Research prompt for literature review and citation assembly.

TOPIC
Evidence that adult or pediatric obesity alters microbial fermentation capacity / SCFA production or fecal SCFA profiles, including studies with null, small, or contradictory effects. Focus on whether obesity status modifies SCFA responses to carbohydrate substrates.

NEEDED OUTPUT
1. Balanced evidence map: supportive, null, and mixed findings.
2. Separate observational fecal SCFA associations from interventional or ex vivo fermentation studies.
3. Sample-size and measurement caveats typical of this literature.
4. Citations that discuss composition vs function discordance (community structure can differ when metabolite output does not, and vice versa).
5. Claim language templates for reporting nonsignificant obesity-group SCFA differences WITHOUT claiming preserved/equivalent capacity.

CONSTRAINTS
Our study: n=8 per group; six SCFA difference-in-change contrasts nonsignificant; equivalence not tested. Do not supply citations framed as proving “no difference = preserved capacity.”
Return: evidence map table + recommended Discussion wording caveats + citations with DOI/PMID.
```

---

### Methods

#### Prompt 06 — Human fecal collection, banking, and ex vivo inoculum preparation standards

**Manuscript target:** Methods 2.1–2.3  
**Citation job:** Best-practice and precedent citations for sample handling.

```text
Research prompt for literature review and citation assembly.

TOPIC
Best practices and published precedents for human stool collection for microbiome + metabolomics studies, including time-to-freeze, storage temperature, freeze–thaw effects on community composition and fermentative activity, and preparation of fecal inocula for anaerobic batch culture.

NEEDED OUTPUT
1. Leading citations / SOPs / reviews on stool biobanking effects relevant to fermentation assays.
2. Citations justifying reporting of antibiotics, diet restrictions, and recent medication exposures.
3. Precedents for splitting stool into paired aliquots (technical biological split) and analyzing nested technical/culture replicates.
4. Recommended reporting items for Methods: atmosphere, temperature, pH control or monitoring, agitation, medium composition, dose of test carbohydrate, inoculum % (w/v), termination method.
5. Note common pitfalls that inflate degrees of freedom when aliquots or wells are treated as independent participants.

CONSTRAINTS
Adolescent focus if literature exists; otherwise state adult precedent with pediatric applicability note. Include DOI/PMID.
Return: Methods citation checklist (item → citation) + 5–8 reporting-standard citations.
```

#### Prompt 07 — Targeted SCFA quantification by PFBBr / GC-MS and interpretive reporting

**Manuscript target:** Methods 2.4; Results/Discussion language  
**Citation job:** Analytical method citations; LOD/LOQ/QC norms; unit and “accumulation” reporting.

```text
Research prompt for literature review and citation assembly.

TOPIC
Quantitative analysis of short-chain fatty acids in biological matrices by pentafluorobenzyl bromide (PFBBr) derivatization with GC-MS under negative chemical ionization (or closely related GC-MS SCFA methods), including calibration, internal standards, QC, LOD/LOQ, and batch handling.

NEEDED OUTPUT
1. Canonical analytical chemistry citations for PFBBr-SCFA workflows and acceptable alternatives commonly used in microbiome studies.
2. Citations on interpreting culture supernatant vs fecal SCFA measurements.
3. Reporting standards for units, dilution correction, and below-LOQ handling.
4. Method papers that motivate reporting succinate and 5-aminovalerate as related fermentation metabolites when included.
5. Language guidance: concentration / accumulation / net output vs production.

CONSTRAINTS
We need facility-confirmable methods citations, not product ads. Prefer peer-reviewed analytical methods over vendor notes.
Return: Methods citation block + 6 reporting best-practice bullets with citations.
```

#### Prompt 08 — 16S V3–V4 processing, genus-level inference limits, and absolute/qPCR-scaled abundance

**Manuscript target:** Methods 2.5–2.6; Discussion 4.6, 4.8  
**Citation job:** DADA2/phyloseq norms; why genus ≠ species function; relative vs absolute abundance.

```text
Research prompt for literature review and citation assembly.

TOPIC
Standard practices and leading citations for:
1) 16S rRNA amplicon sequencing (V3–V4), DADA2 ASV inference, taxonomy assignment, contaminant/depth filters, genus aggregation
2) Limitations of genus-level 16S for attributing metabolic function
3) Total bacterial 16S qPCR used to scale relative abundances to estimated absolute / “qPCR-scaled” abundances
4) Why relative-abundance increases can occur without absolute expansion

NEEDED OUTPUT
1. Canonical methods citations (DADA2, phyloseq, Silva/GTDB taxonomy caveats as relevant).
2. Key papers on compositional data bias and absolute abundance estimation strategies.
3. Citations justifying terminology: “qPCR-scaled genus abundance” / “estimated 16S target equivalents,” and against calling this species-specific cell counts without qualification.
4. Precedents for reporting relative and qPCR-scaled differential abundance as complementary estimands rather than method-voting.

CONSTRAINTS
Our differential abundance used ANCOM-BC2 (relative), MaAsLin3 (qPCR-scaled), and ALDEx3 (scale uncertainty). Prefer citations that justify triangulation without claiming concordance = truth.
Return: Methods citation checklist + Discussion caveats for 16S/qPCR-scaled inference.
```

#### Prompt 09 — Differential abundance frameworks: ANCOM-BC2, MaAsLin, ALDEx triangulation

**Manuscript target:** Methods 2.8; Results 3.7; Discussion 4.6  
**Citation job:** Method papers + interpretive norms for exploratory q thresholds.

```text
Research prompt for literature review and citation assembly.

TOPIC
Statistical frameworks for microbiome differential abundance:
- ANCOM-BC / ANCOM-BC2
- MaAsLin2 / MaAsLin3
- ALDEx2 / ALDEx3
Focus on what each estimand assumes, how compositionality and scale uncertainty are handled, and why methods often disagree.

NEEDED OUTPUT
1. Primary methods papers for each framework (cite exact intended versions when possible).
2. Benchmark / comparison papers showing incomplete concordance across DA tools.
3. Guidance citations on exploratory false-discovery thresholds (e.g., q<0.1) vs confirmatory thresholds.
4. Recommended reporting language when one genus is significant on relative abundance but not on qPCR-scaled absolute-scale models.
5. Citations warning against “method counting” as a significance test.

CONSTRAINTS
Map recommendations to an exploratory head-to-head carbohydrate contrast at a single timepoint with nested repeated measures. Include DOI/PMID.
Return: short methods-rationale paragraph + citation table + 5 interpretive rules.
```

#### Prompt 10 — Nested mixed models, ICC, and donor-aware PERMANOVA for paired fermentation designs

**Manuscript target:** Methods 2.8; Results 3.5–3.6; Discussion strengths  
**Citation job:** Justify donor→aliquot→well hierarchy, ICC as repeatability, and PERMANOVA/betadisper reporting.

```text
Research prompt for literature review and citation assembly.

TOPIC
Statistical citations and precedents for analyzing nested ex vivo fermentation experiments where:
- Donor is the independent biological unit
- Paired aliquots are nested within donor
- Culture wells are nested within aliquot × condition and link repeated timepoints
Outcomes include longitudinal metabolite concentrations analyzed with linear mixed models, plus community Bray–Curtis distances analyzed with PERMANOVA and multivariate dispersion tests (betadisper).

NEEDED OUTPUT
1. Citations / textbooks / papers for nested random effects and estimated marginal means / Tukey-adjusted contrasts.
2. ICC citations for quantifying interpersonal repeatability of biomarker responses across conditions.
3. vegan/PERMANOVA and Anderson dispersion literature; warnings about inferring term significance from omnibus models alone; importance of term-specific R².
4. Precedents for averaging nested technical replicates before person-level microbiome inference.
5. Bootstrap / permutation sensitivity analysis precedents when residual non-normality or heteroscedasticity is present.

CONSTRAINTS
Prefer biostatistics + microbiome ecology sources usable in a Methods paragraph. Include DOI/PMID or standard textbook refs.
Return: Methods citation block structured as LMM / ICC / PERMANOVA / diagnostics.
```

---

### Discussion

#### Prompt 11 — Interpersonal heterogeneity in microbial fermentation as the conceptual center

**Manuscript target:** Discussion 4.5 (centerpiece); Results 3.5 framing  
**Citation job:** Prior work on donor-level response phenotypes and personalization limits.

```text
Research prompt for literature review and citation assembly.

TOPIC
Interpersonal (donor-level) heterogeneity in gut microbial fermentation of identical carbohydrate substrates, including evidence that individual ranking of high-to-low SCFA response is relatively stable across related substrates (generalized responsiveness) versus substrate-specific responder phenotypes.

NEEDED OUTPUT
1. Leading human studies documenting large donor effects on SCFA output in vitro/ex vivo.
2. Evidence for cross-substrate correlation / generalized fermentative responsiveness.
3. Evidence (if any) for stable substrate-specific phenotypes, and how such phenotypes should be validated.
4. Citations cautioning against median-split “responder/non-responder” labels as clinical phenotypes.
5. Claim set supporting future continuous-phenotype framing without claiming our study validates a classifier or dietary predictor.

OUR RESULT TO CONTEXTUALIZE (do not overclaim)
Donor ICCs ~0.63–0.84; RDC–SDC Spearman ρ ~0.83–0.96 for acetate/propionate/butyrate change. A/B median-split agreement only 7–12/16 donors.

CONSTRAINTS
No precision-nutrition diagnostic claims. Include DOI/PMID.
Return: Discussion-ready prior-claim bullets + citation table + competing-interpretation notes.
```

#### Prompt 12 — Background fermentation in no-substrate controls and non-monotonic SCFA responses

**Manuscript target:** Discussion 4.3  
**Citation job:** Explain why no-added-carbohydrate wells still accumulate SCFAs and why some active-vs-control contrasts can be null or negative.

```text
Research prompt for literature review and citation assembly.

TOPIC
In fecal batch-fermentation assays, SCFAs often accumulate even without added test carbohydrate because of residual fermentable material in medium and inoculum. Also compile explanations for analyte-specific or negative contrasts (e.g., lower propionate change under a test carbohydrate vs control), including pH effects, substrate competition, cross-feeding shifts, net consumption/conversion of SCFAs, and community adaptation over 24–48 h.

NEEDED OUTPUT
1. Methods and primary papers documenting nontrivial control-well SCFA accumulation.
2. Mechanistic/interpretive citations for cross-feeding and SCFA utilization.
3. Examples where “added carbohydrate” did not increase all SCFAs uniformly.
4. Claim language for discussing our analyte-specific adjusted contrasts without forcing a more-substrate-more-SCFA narrative.

OUR RESULT TO CONTEXTUALIZE
All conditions including no-added-carbohydrate accumulated acetate, propionate, butyrate over 48 h. Many active-vs-control change contrasts were nonsignificant; RDC propionate change was lower than no-added-carbohydrate in adjusted contrasts. SDC vs RDC differences were nonsignificant for all three primary SCFAs.

CONSTRAINTS
Treat mechanistic explanations as hypotheses unless directly demonstrated in cited papers. Include DOI/PMID.
Return: Discussion paragraph scaffold + citation-backed explanatory hypotheses list.
```

#### Prompt 13 — Composition–function discordance: community structure vs SCFA output

**Manuscript target:** Discussion 4.4–4.5; Results 3.6 framing  
**Citation job:** Why Bray–Curtis can show group/time effects while carbohydrate metabolites show different patterns.

```text
Research prompt for literature review and citation assembly.

TOPIC
Literature on discordance (or concordance) between microbiome community structure metrics (alpha/beta diversity, PERMANOVA) and functional metabolic outputs (SCFA profiles) under dietary or in vitro carbohydrate interventions.

NEEDED OUTPUT
1. Key papers where community shifted without SCFA differences, or SCFAs shifted without global community PERMANOVA signals.
2. Reviews arguing that functional redundancy can stabilize metabolite output despite compositional change.
3. Citations supporting reporting term-specific PERMANOVA R² separately for time, group, substrate, and interactions.
4. Language for stating a modest obesity-group community association alongside nonsignificant SCFA group interactions.

OUR RESULT TO CONTEXTUALIZE
Donor-aggregated PERMANOVA: strong time effect (R²≈0.24); modest group effect (R²≈0.026 overall; ≈0.055 at 48 h); nonsignificant carbohydrate and group×carbohydrate terms at claim-facing models. SCFA obesity difference-in-change contrasts nonsignificant.

CONSTRAINTS
Do not cite sources to claim carbohydrate drove 48 h beta diversity. Include DOI/PMID.
Return: Discussion claim bullets + citation table.
```

#### Prompt 14 — *Fusicatenibacter*: ecology, fermentation associations, and evidence quality

**Manuscript target:** Discussion 4.6; Results 3.7–3.8  
**Citation job:** What is known; what is not; avoid mechanistic overreach.

```text
Research prompt for literature review and citation assembly.

TOPIC
Genus Fusicatenibacter (Firmicutes / Lachnospiraceae-associated lineage as classified in current taxonomy): taxonomic history, typical human gut prevalence, associations with diet or SCFA-related phenotypes, any evidence of carbohydrate utilization or butyrate production at strain level, and quality/limits of that evidence.

NEEDED OUTPUT
1. Foundational taxonomy / description papers.
2. Human observational associations (IBD, metabolic disease, diet interventions)—balanced, noting inconsistency.
3. Any culture-based, genomic, or metabolomic evidence regarding butyrate production or starch/carbohydrate utilization.
4. Explicit statement of whether current literature supports calling Fusicatenibacter a confirmed butyrate producer on our substrates.
5. Recommended cautious Discussion wording for an exploratory relative-abundance SDC>RDC signal that lacks significant qPCR-scaled expansion and lacks a significant donor-level taxon–butyrate association.

OUR RESULT TO CONTEXTUALIZE
Exploratory ANCOM-BC2: relative Fusicatenibacter higher under SDC than RDC at 48 h (q<0.1 exploratory). MaAsLin3/ALDEx3 qPCR-scaled SDC×48 h nonsignificant after multiplicity. Donor-level SDC-attributable net-change Spearman ρ=0.477, permutation p=0.088, bootstrap CI crosses zero (N=14). No isolate/flux experiments.

CONSTRAINTS
Do not assemble a citation chain that converts an exploratory relative marker into a mechanism. Prefer primary data over secondary blog/review hype. Include DOI/PMID.
Return: evidence grades (strong/moderate/weak/absent) for ecology, butyrate production, starch utilization, human metabolic association + Discussion-safe claim bullets.
```

#### Prompt 15 — Translating rodent carbohydrate–microbiome diet studies to human ex vivo systems

**Manuscript target:** Discussion 4.2; Introduction ¶3–4  
**Citation job:** Bridge Plaza-Díaz 2022 and related animal work without false equivalence.

```text
Research prompt for literature review and citation assembly.

TOPIC
Translation challenges from rodent diet-replacement studies of rapid vs slow/resistant carbohydrates (microbiome + fecal SCFA + host metabolic endpoints) to human fecal ex vivo fermentation assays that bypass host digestion.

ANCHOR PAPER
Plaza-Díaz et al., Frontiers in Nutrition 2022, DOI 10.3389/fnut.2022.992682 (and closely related slow-digesting / resistant carbohydrate rodent studies).

NEEDED OUTPUT
1. Summarize what the anchor paper and nearest related animal studies claimed (microbiome pathways, fecal acetate/propionate, host outcomes).
2. Systematic list of design differences that block direct claim transfer (digestion, dose, diet matrix, colonization, host metabolism, fecal vs culture SCFAs, age/species).
3. Citations from translational microbiome methods literature on when animal dietary microbiome findings should generate human hypotheses vs when they should not be framed as replication targets.
4. Discussion language for “biological rationale / hypothesis generator,” not “attempted replication.”

CONSTRAINTS
Our study directly exposes adolescent fecal communities to RDC/SDC and finds nonsignificant SDC–RDC SCFA change differences. Do not retrieve literature to reframe our null head-to-head as a successful translation of rodent superiority claims.
Return: side-by-side design comparison table + Discussion paragraph scaffold + citations.
```

#### Prompt 16 — Limits of concurrent microbiome–metabolite association and required mechanistic follow-up

**Manuscript target:** Discussion 4.6–4.8; Conclusion  
**Citation job:** Causal/predictive standards for next experiments.

```text
Research prompt for literature review and citation assembly.

TOPIC
Standards of evidence for linking specific taxa to SCFA responses in human microbiome studies:
- concurrent correlation vs predictive baseline models
- relative vs absolute abundance
- need for isolate growth, substrate disappearance, isotope tracing, metatranscriptomics, and co-culture
Compile authoritative methods/review citations that define what is required before claiming a taxon utilizes a substrate or produces an observed metabolite increase.

NEEDED OUTPUT
1. Review papers on microbiome causality frameworks (Koch/ Bradford-Hill adapted; multi-omics triangulation).
2. Exemplar studies that successfully (or carefully) linked taxa to SCFA via mechanism—not just correlation.
3. Checklist of follow-up experiments appropriate after an exploratory relative-abundance hit (our Fusicatenibacter case).
4. Citations cautioning against predictive claims from concurrent same-experiment features without held-out validation.

CONSTRAINTS
Our manuscript must nominate candidates for follow-up, not claim mechanism or prediction. Include DOI/PMID.
Return: “claims allowed now” vs “claims requiring new data” table + prioritized follow-up citation list.
```

#### Prompt 17 — Adolescent / pediatric gut microbiome fermentation literature specifically

**Manuscript target:** Introduction ¶1 & ¶4; Discussion strengths/limitations  
**Citation job:** Age-specific citations to justify adolescent sampling and state pediatric evidence gaps.

```text
Research prompt for literature review and citation assembly.

TOPIC
Pediatric and adolescent gut microbiome literature relevant to obesity and carbohydrate fermentation, with emphasis on:
1) developmental differences vs adult microbiota
2) any ex vivo / in vitro fermentation studies using child or adolescent stool
3) gaps that make adult-only literature insufficient for adolescent claims

NEEDED OUTPUT
1. Key developmental microbiome reviews (age ~10–18 if possible).
2. Pediatric obesity microbiome studies that measure function (SCFAs, fermentation) rather than only relative taxa.
3. Explicit gap statement suitable for Introduction translational paragraph.
4. Ethics/consent/assent reporting precedents for adolescent stool research (brief).

CONSTRAINTS
If pediatric fermentation studies are scarce, document that scarcity with search evidence rather than substituting adult findings silently. Include DOI/PMID.
Return: age-stratified citation table + gap statement + 5 Introduction-safe claims.
```

#### Prompt 18 — Competitive prior art: human trials and fermentation studies of slow/rapid digestible carbohydrates

**Manuscript target:** Introduction ¶3; Discussion 4.2–4.3  
**Citation job:** Position the paper among nearest human evidence without overstating novelty incorrectly.

```text
Research prompt for literature review and citation assembly.

TOPIC
Identify the closest prior human research on rapid-digesting vs slow-digesting / slowly digestible / resistant starch carbohydrate formulations and microbial outcomes, including:
- randomized or controlled feeding studies measuring microbiome and/or SCFAs
- human fecal in vitro fermentation of SDS/RS vs rapidly digestible starch
- commercial multi-component slow-digesting blends if published independently of marketing materials

NEEDED OUTPUT
1. Shortlist of 10–15 nearest neighbors with design sketches (population, substrate, endpoints, key result).
2. Novelty/positioning statement templates for our adolescent paired ex vivo RDC/SDC design with obesity stratification and donor-repeatability quantification.
3. Claims others made that our claim lock forbids us from copying.
4. Gaps our study fills (and gaps it does not).

CONSTRAINTS
Prefer peer-reviewed human data. Flag industry-affiliated studies transparently. Include DOI/PMID.
Return: competitive landscape table + 1-paragraph positioning text + citation list.
```

---

## Recommended run order

| Order | Prompt | Why first |
|---|---|---|
| 1 | 02, 04 | Locks correct SCFA / assay language for all later drafting |
| 2 | 03, 15, 18 | Positions against rodent anchor + human prior art |
| 3 | 01, 05, 17 | Completes Introduction funnel |
| 4 | 06–10 | Builds Methods citation spine |
| 5 | 11–14, 16 | Writes Discussion centerpiece and *Fusicatenibacter* caution |

## Downstream synthesis format (after LLM runs)

For manuscript drafting, condense each prompt’s output into:

```text
SECTION:
PRIOR CLAIM:
CITATIONS (max 3 lead + 2 backup):
CONFLICTS / NUANCE:
MANUSCRIPT-SAFE SENTENCE:
CLAIM-LOCK CHECK (supported / exploratory / unsupported if used beyond evidence):
```

## Explicit non-goals for this citation package

These prompts are **not** meant to justify:

- SDC superiority for acetate, propionate, or butyrate
- obesity-group SCFA equivalence or “preserved capacity”
- carbohydrate-driven 48-hour Bray–Curtis restructuring
- absolute *Fusicatenibacter* expansion
- *Fusicatenibacter* as the mechanistic butyrate source
- validated responder phenotypes or dietary-response prediction
- clinical, glycemic, or obesity-treatment benefit of RDC/SDC
