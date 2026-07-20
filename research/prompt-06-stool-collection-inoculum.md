# Prompt 06 — Stool collection, biobanking, and fecal inoculum preparation for anaerobic batch culture

**Manuscript target:** Methods 2.1–2.3 (donor recruitment/screening, stool collection & storage, fermentation setup) · **Citation job:** supply verified methodological precedent and reporting-standard citations for stool handling, freeze–thaw, paired-aliquot splitting, batch-culture parameter reporting, and pseudoreplication avoidance · **Date:** 2026-07-19

All PMIDs/DOIs below were retrieved from PubMed (via the bio-research PubMed MCP) or Consensus and cross-checked; Hurlbert 1984 was verified via the ESA/Wiley DOI landing page. No citation in this report is reconstructed from memory.

---

## Prioritized citation table

| # | Title | Year | Journal | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 1 | An in vitro batch fermentation protocol for studying the contribution of food to gut microbiota composition and functionality | 2021 | Nature Protocols | [10.1038/s41596-021-00537-x](https://doi.org/10.1038/s41596-021-00537-x) · PMID 34089022 | **Lead Methods protocol citation.** Canonical published batch-fermentation SOP specifying fermentation medium (oligotrophic, peptone-based), fecal inoculum concentration, anaerobic handling, incubation time, and termination. Anchors 2.3. |
| 2 | INFOGEST static in vitro simulation of gastrointestinal food digestion | 2019 | Nature Protocols | [10.1038/s41596-018-0119-1](https://doi.org/10.1038/s41596-018-0119-1) · PMID 30886367 | International consensus upstream digestion protocol that #1 extends. Cite if RDC/SDC substrates were pre-digested before fermentation; also the standard reference for why digestion conditions must be reported explicitly. |
| 3 | Evaluation of an optimal preparation of human standardized fecal inocula for in vitro fermentation studies | 2015 | J Microbiol Methods | [10.1016/j.mimet.2015.07.019](https://doi.org/10.1016/j.mimet.2015.07.019) · PMID 26222994 | **Core freeze-vs-fresh inoculum citation.** Directly compared four fecal preservation routes against fresh feces for SCFA output and community composition; glycerol + dialysate + −80 °C best preserved fermentative activity. Bacteroidetes diversity declined under all frozen treatments — the key caveat to state. |
| 4 | Methods for Improving Human Gut Microbiome Data by Reducing Variability through Sample Processing and Storage of Stool | 2015 | PLoS ONE | [10.1371/journal.pone.0134802](https://doi.org/10.1371/journal.pone.0134802) · PMID 26252519 | **Core time-to-freeze + homogenization + freeze–thaw citation.** Empirical basis for freezing within ~15 min, limiting frost-free freezer storage to <3 days, homogenizing whole stool before subsampling, and the finding that freeze–thaw altered taxon abundance only beyond ~4 cycles. |
| 5 | Storage conditions of intestinal microbiota matter in metagenomic analysis | 2012 | BMC Microbiology | [10.1186/1471-2180-12-158](https://doi.org/10.1186/1471-2180-12-158) · PMID 22846661 | Widely cited primary demonstration that storage temperature/duration before DNA extraction shifts metagenomic profiles; standard companion citation to #4. |
| 6 | Sample pre-treatment procedures for the omics analysis of human gut microbiota: turning points, tips and tricks for gene sequencing and metabolomics | 2020 | J Pharm Biomed Anal | [10.1016/j.jpba.2020.113592](https://doi.org/10.1016/j.jpba.2020.113592) · PMID 32947167 | **Best single review** covering the combined 16S + metabolomics preanalytical chain: collection, transport, storage, pretreatment, and the confounders (antibiotics, diet, disease) that must be recorded. Justifies both 2.1 screening items and 2.2 handling items. |
| 7 | Short-Chain Fatty Acid Production by Gut Microbiota from Children with Obesity Differs According to Prebiotic Choice and Bacterial Community Composition | 2020 | mBio | [10.1128/mBio.00914-20](https://doi.org/10.1128/mBio.00914-20) · PMID 32788375 | **Closest design precedent, and adolescent-specific.** Ex vivo fecal fermentation of carbohydrate substrates using stool from 17 children/adolescents with obesity (ages 10–18), SCFA endpoints + 16S, with explicit donor-dependent (interdonor) variation in SCFA production. Establishes both the pediatric applicability and the donor-as-unit framing. |
| 8 | Quantifying spatiotemporal variability and noise in absolute microbiota abundances using replicate sampling (DIVERS) | 2019 | Nature Methods | [10.1038/s41592-019-0467-y](https://doi.org/10.1038/s41592-019-0467-y) · PMID 31308552 | **Lead citation for nested replicate variance decomposition.** Formal method partitioning technical noise vs. spatial (within-stool) sampling variability vs. biological/temporal variability; found technical noise dominated abundance variability for nearly half of gut taxa. Directly supports treating aliquots/wells as nested, not independent. |
| 9 | Fecal Aliquot Straw Technique (FAST) allows for easy and reproducible subsampling | 2018 | Microbiome | [10.1186/s40168-018-0458-8](https://doi.org/10.1186/s40168-018-0458-8) · PMID 29776435 | **Lead precedent for paired-aliquot splitting.** Explicitly split one stool into replicate straw aliquots and showed 16S community reproducibility between paired aliquots (discordant OTUs <0.3% abundance), while preserving viable microbes for downstream culture — exactly the technical-biological split logic. |
| 10 | What exactly is 'N' in cell culture and animal experiments? | 2018 | PLoS Biology | [10.1371/journal.pbio.2005282](https://doi.org/10.1371/journal.pbio.2005282) · PMID 29617358 | **Lead pseudoreplication citation.** Distinguishes biological units, experimental units, and observational units; documents 46% pseudoreplication rate in a random published sample; gives explicit *ex vivo* and *in vitro* design examples. Perfect fit for "donor = experimental unit, well = observational unit." |
| 11 | Pseudoreplication and the design of ecological field experiments | 1984 | Ecological Monographs | [10.2307/1942661](https://doi.org/10.2307/1942661) *(no PMID; not indexed in PubMed — verified via ESA/Wiley DOI landing page)* | Seminal definition of pseudoreplication and of inflated degrees of freedom from non-independent replicates. Optional historical anchor alongside #10. |
| 12 | Reporting guidelines for human microbiome research: the STORMS checklist | 2021 | Nature Medicine | [10.1038/s41591-021-01552-x](https://doi.org/10.1038/s41591-021-01552-x) · PMID 34789871 | **Primary reporting standard.** 17-item checklist across six manuscript sections, including specimen collection/handling, laboratory, bioinformatics, and statistical reporting items. Cite as the framework the Methods section follows. |
| 13 | Minimum information about a marker gene sequence (MIMARKS) and minimum information about any (x) sequence (MIxS) specifications | 2011 | Nature Biotechnology | [10.1038/nbt.1823](https://doi.org/10.1038/nbt.1823) · PMID 21552244 | Reporting standard for the 16S V3–V4 amplicon metadata (primers, region, environmental package, sample provenance). Complements #12 for the sequencing half of Methods 2.2. |
| 14 | Assessment of variation in microbial community amplicon sequencing by the Microbiome Quality Control (MBQC) project consortium | 2017 | Nature Biotechnology | [10.1038/nbt.3981](https://doi.org/10.1038/nbt.3981) · PMID 28967885 | Multi-laboratory quantification showing variability depended most on specimen type/origin, then DNA extraction, sample-handling environment, and bioinformatics. Justifies reporting handling and extraction details as first-class Methods items. |
| 15 | Population-level analysis of gut microbiome variation | 2016 | Science | [10.1126/science.aad3503](https://doi.org/10.1126/science.aad3503) · PMID 27126039 | **Lead justification for recording covariates.** 69 covariates replicably associated with microbiota composition; stool consistency had the largest effect size and medication explained the largest total variance. Explicit call to include host covariates in study design. |
| 16 | Impact of commonly used drugs on the composition and metabolic function of the gut microbiota | 2020 | Nature Communications | [10.1038/s41467-019-14177-z](https://doi.org/10.1038/s41467-019-14177-z) · PMID 31953381 | Justifies recording *non-antibiotic* medication exposures (PPIs, metformin, laxatives, antibiotics strongest), and explicitly demonstrates the need to correct for multiple drug use. |
| 17 | Recovery of gut microbiota of healthy adults following antibiotic exposure | 2018 | Nature Microbiology | [10.1038/s41564-018-0257-9](https://doi.org/10.1038/s41564-018-0257-9) · PMID 30349083 | **Empirical basis for an antibiotic washout window.** Composition recovered to near-baseline within ~1.5 months after a 4-day broad-spectrum intervention, but 9 common species remained undetectable at day 180; butyrate producers were depleted early. Supports (and bounds) whatever exclusion window the protocol used. |
| 18 | Diet rapidly and reproducibly alters the human gut microbiome | 2013 | Nature | [10.1038/nature12820](https://doi.org/10.1038/nature12820) · PMID 24336217 | Justifies recording/standardizing recent diet: short-term macronutrient change altered community structure within days and shifted the carbohydrate-vs-protein fermentation balance. |
| 19 | The effect of in vitro simulated colonic pH gradients on microbial activity and metabolite production using common prebiotics as substrates | 2024 | BMC Microbiology | [10.1186/s12866-024-03235-2](https://doi.org/10.1186/s12866-024-03235-2) · PMID 38468200 | **Lead pH-reporting citation.** Shows pH regime significantly changes in vitro community structure and SCFA output in a donor- and substrate-dependent manner — and that donor and substrate effects were larger than pH. Justifies reporting pH control/monitoring *and* supports donor as the dominant variance source. |
| 20 | GC-MS method development and validation for the determination of short chain fatty acids in human feces | 2026 | J Pharm Biomed Anal | [10.1016/j.jpba.2026.117488](https://doi.org/10.1016/j.jpba.2026.117488) · PMID 41967141 | Recent validated SCFA quantification method with explicit stability data: SCFAs stable in acidified fecal samples ~10 days without cold chain; −80 °C optimal for long-term, 4 °C for short-term. Supports the termination/acidification and analyte-storage sentences. |
| 21 | Impact of fructooligosaccharides on gut microbiota composition and metabolite production: implications for childhood obesity | 2025 | PeerJ | [10.7717/peerj.19894](https://doi.org/10.7717/peerj.19894) · PMID 40895066 | Second pediatric precedent: ex vivo fermentation of stool from 39 children aged 6–15 with/without a carbohydrate substrate, SCFA and gas endpoints. Useful for "pediatric donors have been used in this assay format." |
| 22 | Evaluating gut microbiota profiles from archived fecal samples | 2018 | BMC Gastroenterology | [10.1186/s12876-018-0896-6](https://doi.org/10.1186/s12876-018-0896-6) · PMID 30409123 | 16S V3–V4-specific storage/archival assessment; found short-term room-temperature stability but significantly altered diversity in long-term archived samples. Supports stating archive duration explicitly. |

---

## Prior-claim bullets (each mapped to citation numbers)

**Stool biobanking effects relevant to fermentation assays**

- Time-to-freeze is a documented, measurable source of microbiome variability; freezing within ~15 minutes of defecation and limiting frost-free freezer storage to <3 days have been recommended on empirical grounds [4]. Room-temperature stability over shorter intervals has also been reported [22], so the literature is *directionally consistent but not numerically unanimous* — see Conflicts.
- Storage temperature and storage duration before nucleic-acid extraction alter recovered community profiles [5, 4, 22].
- Freeze–thaw cycling is a real but comparatively *weak* perturbation of DNA-based composition relative to time-to-freeze and subsampling: an effect on taxon abundance was detected only beyond roughly four cycles [4].
- Freeze–thaw effects on **fermentative activity** are a separate question from effects on **composition**, and the fermentation-specific evidence is the load-bearing citation: cryopreserved fecal inocula retained the qualitative SCFA profile (acetate > butyrate/propionate) but the *quantitative* SCFA output differed by preservation method, and Bacteroidetes diversity declined under all frozen treatments [3].
- Whole-stool homogenization prior to subsampling substantially reduces intra-sample variability, because a single stool is spatially heterogeneous [4, 8].
- Preanalytical handling (specimen origin, extraction, handling environment) contributes more variance to amplicon profiles than bioinformatics choices [14], and this is true across the combined sequencing + metabolomics workflow [6].

**Justifying reporting of antibiotics, diet, and recent medication exposures**

- Medication use explains the largest total share of microbiota compositional variance among measured covariates in population cohorts, and stool consistency has the single largest individual effect size; the authors explicitly urge inclusion of host covariates in study design [15].
- Non-antibiotic drugs (proton-pump inhibitors, metformin, laxatives) show among the strongest microbiome associations, and multiple concurrent drug use must be adjusted for [16].
- Antibiotic exposure depletes butyrate producers and *Bifidobacterium* early, with near-baseline compositional recovery by ~1.5 months but incomplete species-level recovery at 180 days [17]. This bounds any washout window claim in either direction.
- Short-term diet change alters community structure and the carbohydrate-vs-protein fermentation balance within days [18], which is directly relevant to a carbohydrate-fermentation endpoint.
- Preanalytical reviews list antibiotics, diet, age, and disease among the confounders that should be captured at collection [6].

**Paired aliquots and nested technical/culture replicates**

- Splitting one stool into replicate aliquots and demonstrating between-aliquot concordance is an established, published technique; paired aliquots from the same stool yielded highly reproducible 16S communities, with discordant OTUs confined to taxa <0.3% relative abundance, while preserving microbial viability for downstream culture [9].
- Formal variance decomposition across replicate samples — separating technical noise, within-specimen spatial sampling variability, and biological variability — is an established analysis for exactly this nesting structure, and it found technical noise dominant for nearly half of detected gut taxa [8].
- Within-donor, multi-substrate ex vivo fermentation with donor-level interpretation of SCFA differences is the established design in the adolescent-obesity context [7].

**Batch-culture Methods reporting items**

- Published batch-fermentation protocols specify medium composition, inoculum concentration, anaerobic handling, incubation duration, and sampling/termination as required, reportable parameters [1]; upstream digestion conditions are separately specified when substrates are pre-digested [2].
- pH is a first-order determinant of in vitro SCFA output and community structure, and should therefore be controlled or at minimum monitored and reported — though donor and substrate effects exceeded pH effects in the study that tested this directly [19].
- SCFA analyte handling after termination (acidification, storage temperature, hold time) has published stability bounds that should be reported alongside the quantification method [20].

**Pseudoreplication**

- Treating wells, aliquots, or technical replicates as independent participants inflates the effective sample size and the degrees of freedom, producing false positives; the correct framing distinguishes the *biological unit* (donor), the *experimental unit* (the donor × carbohydrate incubation to which the treatment was applied), and the *observational unit* (the well/aliquot measured) [10, 11].
- Empirically, pseudoreplication is common: only 22% of a random sample of published animal studies replicated the correct entity–intervention pair; 46% were pseudoreplicated [10].
- STORMS requires explicit reporting of the statistical model and the unit of analysis [12].

---

## Conflicts / nuance

1. **Time-to-freeze evidence is not unanimous.** Gorzelak et al. recommend freezing within ~15 minutes [4], whereas archival/FIT-based work reports "room temperature stability of the gut microbiome" over the shorter handling intervals tested [22], and a validated SCFA method reports stability in *acidified* feces for up to 10 days without cold chain [20]. These are not contradictory so much as endpoint-specific: DNA-based taxon abundance, archived-sample diversity, and chemically stabilized SCFA analytes have different tolerances. Report the actual observed interval; do not cite a 15-minute rule as universal.

2. **Composition-preserving ≠ activity-preserving.** Freeze–thaw is mild for DNA-based composition below ~4 cycles [4] but demonstrably alters quantitative SCFA output from a cryopreserved inoculum [3]. For a study whose primary endpoints are SCFA accumulation, [3] is the governing citation, not [4]. Bacteroidetes loss on freezing [3] is the specific, named caveat.

3. **Frozen inocula remain a deviation from fresh.** Even the best-performing preservation route in [3] was described as "high similarities to" — not equivalence with — fresh feces. Any manuscript sentence should say "comparable to" or "an accepted alternative to," never "equivalent to."

4. **Donor effect generally exceeds assay-condition effects.** Both the pH-gradient study [19] and the adolescent prebiotic study [7] found donor identity to be a dominant or larger source of variation than the manipulated condition. This *supports* donor-aware modeling but is not itself evidence about carbohydrate-type differences.

5. **Adolescent evidence is thin and confined to the fermentation-design literature, not the biobanking literature.** Pediatric/adolescent ex vivo fermentation precedent exists and is strong [7, 21]. The stool storage, freeze–thaw, and inoculum-preservation evidence [3, 4, 5, 14, 22] is essentially all adult-derived. Recommended manuscript stance: cite the adult preanalytical literature as the methodological basis and add an explicit pediatric applicability note — the preanalytical mechanisms (spatial heterogeneity, cold-chain kinetics, freeze-thaw effects on viability) are not known to be age-specific, but they have not been separately validated in adolescents.

6. **Antibiotic washout has no consensus threshold.** [17] shows recovery is neither immediate nor complete: near-baseline at ~1.5 months, still species-incomplete at 180 days. Whatever exclusion window the study used should be stated as a pragmatic protocol decision citing [17], not as a validated sufficiency threshold.

7. **STORMS was written for observational human microbiome studies**, not ex vivo fermentation assays. It is the right framework for the donor-recruitment, specimen-handling, sequencing, and statistical-reporting items, but the fermentation parameters themselves must be reported per the protocol literature [1, 2, 19]. There is no single consensus reporting checklist for anaerobic batch culture — this is a genuine gap, and the honest framing is "reported in accordance with [1] and [12]."

---

## What we must not overclaim (given claim lock)

- **These are methods citations only.** None of #1–22 supports SDC superiority over RDC for acetate, propionate, or butyrate. [7] and [21] establish that carbohydrate substrates are fermented by pediatric donor microbiota and that responses are donor-dependent — they do **not** rank slow- vs rapid-digesting formulations.
- **Do not use [7] to claim obesity-group equivalence or "preserved fermentation capacity."** [7] reports interdonor variation among children with obesity and explicitly notes that fecal SCFA concentration, microbiota SCFA production capacity, and obesity markers did **not** positively correlate. That is a null/decoupling result about correlation structure, not evidence of group equivalence.
- **Do not use [7] or [21] to claim clinical, glycemic, or obesity-treatment benefit.** [7] frames its own findings as an in vitro hypothesis ("suggest the hypothesis that OTC prebiotic supplements may be unequal…"). Preserve that hedging.
- **Do not use [19] to claim carbohydrate-driven community restructuring.** [19] reports that pH and substrate influence community structure in a continuous in vitro colon setup — this is not evidence for 48-hour Bray–Curtis restructuring in our batch design, and [19]'s own finding is that donor effects exceeded substrate effects.
- **Do not use any citation here to support *Fusicatenibacter* claims, absolute-abundance expansion, responder phenotypes, or microbiome-based prediction of dietary response.** No citation in this report addresses those.
- **Do not present frozen-inoculum use as validated-equivalent to fresh.** [3] is a limitation citation as much as a justification citation, and should appear in both Methods and Limitations.
- **Do not cite [12]/[13]/[14] as evidence that our results are reliable.** They are reporting and QC frameworks; adherence to a checklist is a transparency claim, not a validity claim.

---

## Methods citation checklist (item → citation)

Use this as the mapping from each reportable Methods element to its supporting citation. Items marked **(state, no citation needed)** are pure protocol disclosure.

### 2.1 — Donor recruitment, screening, and eligibility

| Methods item to report | Citation |
|---|---|
| Rationale for recording/excluding recent antibiotic exposure; washout window | [17], [15], [16] |
| Rationale for recording non-antibiotic medications (PPIs, metformin, laxatives) and polypharmacy adjustment | [16], [15] |
| Rationale for recording habitual/recent diet and any pre-collection dietary restriction | [18], [6] |
| Rationale for recording stool consistency (largest single covariate effect size) | [15] |
| Rationale for recording age/BMI category and other host covariates | [15], [6] |
| Adolescent donor precedent for this assay format (ages 10–18; ages 6–15) | [7], [21] |
| Ethics/consent, inclusion & exclusion criteria as reportable items | [12] |

### 2.2 — Stool collection, transport, storage, and aliquoting

| Methods item to report | Citation |
|---|---|
| Collection-to-freeze interval (state the actual elapsed time) | [4], [22] |
| Transport temperature and container/preservative (or absence of preservative) | [4], [6], [14] |
| Storage temperature and total storage duration before use | [5], [4], [22] |
| Number of freeze–thaw cycles experienced by each aliquot | [4], [3] |
| Whole-stool homogenization prior to subsampling; rationale (spatial heterogeneity) | [4], [8] |
| Paired/split aliquot scheme — one stool split into technical-biological aliquots | [9], [8] |
| Cryoprotectant (e.g., glycerol) and resuspension buffer for the inoculum aliquot | [3] |
| DNA extraction method and handling environment as a declared variance source | [14], [12] |
| 16S V3–V4 primers, region, and amplicon metadata | [13], [12], [22] |
| Total bacterial 16S qPCR scaling method | [12] *(report per checklist; no dedicated standard identified — see gap note)* |

### 2.3 — Anaerobic batch fermentation

| Methods item to report | Citation |
|---|---|
| Overall batch-fermentation protocol and its provenance | [1] |
| Upstream in vitro digestion of substrates (if performed) | [2], [1] |
| **Atmosphere** — anaerobic gas mix, degassing/pre-reduction, workstation vs. sealed vessel | [1] |
| **Temperature** — incubation temperature | [1] |
| **pH** — controlled or monitored; starting pH; buffering; endpoint pH | [19], [1] |
| **Agitation** — presence/absence, speed, mode | [1] |
| **Medium composition** — basal medium, peptone, reducing agent, minerals, vitamins | [1] |
| **Dose of test carbohydrate** — mass or % w/v of RDC/SDC and of the no-carbohydrate control | [1], [19] |
| **Inoculum % (w/v)** — fecal slurry concentration and dilution scheme | [1], [3] |
| **Fresh vs. cryopreserved inoculum** and its documented effect on SCFA output | [3] |
| **Incubation duration and sampling timepoints** (0 h, 48 h) | [1], [7] |
| **Termination method** — quenching, acidification, centrifugation, snap-freezing | [1], [20] |
| **SCFA quantification** — method, internal standards, LOD/LOQ, post-termination analyte stability | [20] |
| Number of culture replicates per donor × condition, and their nesting | [8], [9], [10] |
| Blank/no-inoculum and no-carbohydrate control wells | [1] |

### 2.4 — Statistical unit of analysis (cross-referenced from Methods 2.3 to Statistics)

| Methods item to report | Citation |
|---|---|
| Explicit declaration: donor = biological/experimental unit; well = observational unit | [10], [11] |
| Donor-aware modeling (random effect for donor / donor-nested replicates) rather than treating wells as independent n | [10], [8] |
| Variance decomposition across nested technical and biological levels | [8] |
| Statement that technical replicates were averaged or modeled, not counted toward n | [10], [11] |
| Overall statistical reporting per checklist | [12] |

**Gap note.** No consensus reporting checklist specific to anaerobic batch fecal fermentation was identified across these searches. The defensible framing is: "Fermentation conditions are reported following the batch-fermentation protocol of Pérez-Burillo et al. [1]; study-level reporting follows STORMS [12] and sequence metadata follows MIxS/MIMARKS [13]."

---

## Reporting-standard citations (7)

These are the citations to gather in a single "reporting standards" sentence or supplementary statement:

1. **STORMS checklist** — Mirzayi et al. 2021, Nat Med — [10.1038/s41591-021-01552-x](https://doi.org/10.1038/s41591-021-01552-x) · PMID 34789871
2. **MIMARKS / MIxS** — Yilmaz et al. 2011, Nat Biotechnol — [10.1038/nbt.1823](https://doi.org/10.1038/nbt.1823) · PMID 21552244
3. **MBQC preanalytical variance** — Sinha et al. 2017, Nat Biotechnol — [10.1038/nbt.3981](https://doi.org/10.1038/nbt.3981) · PMID 28967885
4. **Batch-fermentation protocol** — Pérez-Burillo et al. 2021, Nat Protoc — [10.1038/s41596-021-00537-x](https://doi.org/10.1038/s41596-021-00537-x) · PMID 34089022
5. **INFOGEST digestion consensus** — Brodkorb et al. 2019, Nat Protoc — [10.1038/s41596-018-0119-1](https://doi.org/10.1038/s41596-018-0119-1) · PMID 30886367
6. **Preanalytical omics review (16S + metabolomics)** — Zubeldia-Varela et al. 2020, J Pharm Biomed Anal — [10.1016/j.jpba.2020.113592](https://doi.org/10.1016/j.jpba.2020.113592) · PMID 32947167
7. **Unit-of-analysis / pseudoreplication standard** — Lazic et al. 2018, PLoS Biol — [10.1371/journal.pbio.2005282](https://doi.org/10.1371/journal.pbio.2005282) · PMID 29617358

---

## Downstream synthesis blocks

---

**SECTION:** Methods 2.1 — Donor screening and covariate capture

**PRIOR CLAIM:** Host covariates — especially medication use, recent antibiotic exposure, diet, and stool consistency — are established determinants of fecal microbiota composition, and published guidance calls for capturing them at enrollment.

**CITATIONS (max 3 lead + 2 backup):** Lead — Falony 2016 (PMID 27126039); Vich Vila 2020 (PMID 31953381); Palleja 2018 (PMID 30349083). Backup — David 2014 (PMID 24336217); Zubeldia-Varela 2020 (PMID 32947167).

**CONFLICTS / NUANCE:** No consensus antibiotic washout threshold exists; Palleja shows near-baseline compositional recovery at ~1.5 months but incomplete species-level recovery at 180 days. All four covariate papers are adult cohorts.

**MANUSCRIPT-SAFE SENTENCE:** "Because medication use, recent antibiotic exposure, habitual diet, and stool consistency are established sources of variation in fecal microbiota composition, these exposures were recorded at enrollment and are reported for all donors."

**CLAIM-LOCK CHECK:** Supported. Descriptive covariate-capture rationale only; makes no claim about group differences, carbohydrate effects, or clinical outcome.

---

**SECTION:** Methods 2.2 — Stool collection, storage, and aliquoting

**PRIOR CLAIM:** Time-to-freeze, storage temperature and duration, freeze–thaw cycling, and within-stool spatial heterogeneity are documented, quantifiable sources of variability in fecal microbiome data, and whole-stool homogenization before subsampling reduces intra-sample variability.

**CITATIONS (max 3 lead + 2 backup):** Lead — Gorzelak 2015 (PMID 26252519); Cardona 2012 (PMID 22846661); Sinha 2017 MBQC (PMID 28967885). Backup — Zubeldia-Varela 2020 (PMID 32947167); Rounge 2018 (PMID 30409123).

**CONFLICTS / NUANCE:** Recommended handling intervals differ by endpoint (DNA-based taxon abundance vs. archived diversity vs. chemically stabilized SCFA analytes); report the observed interval rather than invoking a universal rule. All storage evidence is adult-derived — add a pediatric applicability note.

**MANUSCRIPT-SAFE SENTENCE:** "Handling parameters known to affect fecal microbiome and metabolite measurements — interval from defecation to freezing, storage temperature and duration, number of freeze–thaw cycles, and whether stool was homogenized before subsampling — are reported below, since preanalytical handling is a leading contributor of variance in amplicon-based profiling."

**CLAIM-LOCK CHECK:** Supported.

---

**SECTION:** Methods 2.2 — Paired aliquot split

**PRIOR CLAIM:** Splitting a single stool into replicate aliquots for parallel downstream use is an established technique, with published demonstration that paired aliquots from one stool yield highly reproducible 16S communities while preserving microbial viability for culture.

**CITATIONS (max 3 lead + 2 backup):** Lead — Romano 2018 FAST (PMID 29776435); Gorzelak 2015 (PMID 26252519); Ji 2019 DIVERS (PMID 31308552). Backup — Sinha 2017 (PMID 28967885).

**CONFLICTS / NUANCE:** Between-aliquot reproducibility in Romano was demonstrated for 16S composition, not for fermentative output; concordance for SCFA endpoints from split aliquots is not directly established by that citation. Homogenization before splitting is what makes the aliquots comparable [Gorzelak].

**MANUSCRIPT-SAFE SENTENCE:** "Each donor's stool was homogenized and divided into paired aliquots, a subsampling approach that has been shown to yield reproducible community profiles between aliquots of the same specimen while retaining viable microbiota for downstream culture."

**CLAIM-LOCK CHECK:** Supported, with the caveat that "reproducible" refers to community composition; do not extend the word to SCFA output without our own data.

---

**SECTION:** Methods 2.3 — Fecal inoculum preparation (fresh vs. cryopreserved)

**PRIOR CLAIM:** Cryopreserved fecal inocula are an accepted alternative to fresh feces for in vitro fermentation, with glycerol-containing resuspension and −80 °C storage identified as the best-performing route; however, quantitative SCFA output differs by preservation method and Bacteroidetes diversity declines under freezing.

**CITATIONS (max 3 lead + 2 backup):** Lead — Aguirre 2015 (PMID 26222994); Pérez-Burillo 2021 (PMID 34089022). Backup — Gorzelak 2015 (PMID 26252519).

**CONFLICTS / NUANCE:** This is simultaneously a justification and a limitation citation. Aguirre reports "high similarities to" fresh feces, not equivalence. Cite it in both Methods and Limitations.

**MANUSCRIPT-SAFE SENTENCE:** "Fecal inocula were prepared from cryopreserved material, an approach previously evaluated against fresh feces and proposed as an acceptable alternative for in vitro fermentation; we note that preservation method has been reported to affect the quantitative magnitude of short-chain fatty acid production and to reduce recovered Bacteroidetes diversity, and we therefore applied identical handling to every donor and condition."

**CLAIM-LOCK CHECK:** Supported, provided "acceptable alternative" is not upgraded to "equivalent." The identical-handling clause is the honest mitigation and should be retained.

---

**SECTION:** Methods 2.3 — Fermentation conditions reporting

**PRIOR CLAIM:** Published batch-fermentation protocols require explicit reporting of atmosphere, temperature, medium composition, inoculum concentration, substrate dose, incubation duration, and termination method; pH is a first-order determinant of in vitro SCFA output and should be controlled or monitored and reported.

**CITATIONS (max 3 lead + 2 backup):** Lead — Pérez-Burillo 2021 (PMID 34089022); Xie 2024 (PMID 38468200); Brodkorb 2019 (PMID 30886367). Backup — Gkanali 2026 (PMID 41967141).

**CONFLICTS / NUANCE:** No consensus reporting checklist exists specifically for anaerobic batch fecal fermentation. Xie 2024 found donor and substrate effects exceeded pH effects — cite it for the reporting obligation, not as evidence of substrate ranking.

**MANUSCRIPT-SAFE SENTENCE:** "Fermentation conditions — incubation atmosphere, temperature, medium composition, agitation, starting and terminal pH, carbohydrate dose, inoculum concentration (% w/v), incubation duration, and termination procedure — are reported in full, since each has been shown to influence in vitro microbial metabolite output."

**CLAIM-LOCK CHECK:** Supported.

---

**SECTION:** Methods 2.3 / Statistics — Unit of analysis and pseudoreplication

**PRIOR CLAIM:** In ex vivo designs, the donor is the biological and experimental unit and individual culture wells are observational units; treating wells or aliquots as independent participants inflates the effective sample size and degrees of freedom and generates false positives. Nested technical and biological variance can be formally decomposed.

**CITATIONS (max 3 lead + 2 backup):** Lead — Lazic 2018 (PMID 29617358); Ji 2019 (PMID 31308552); Hurlbert 1984 (DOI 10.2307/1942661). Backup — Mirzayi 2021 STORMS (PMID 34789871).

**CONFLICTS / NUANCE:** Lazic's empirical survey is of animal experiments, though the paper explicitly provides in vitro and ex vivo design examples, which is the relevant transfer. Ji's DIVERS uses spike-in absolute abundance quantification, so its variance-partitioning machinery is not directly transplantable to relative-abundance 16S data — cite it for the conceptual nesting argument and for the empirical finding that technical noise is large, not as the specific estimator used.

**MANUSCRIPT-SAFE SENTENCE:** "Because replicate culture wells derived from a single donor's stool are not independent biological replicates, the donor was treated as the experimental unit throughout; technical and culture replicates were modeled as nested within donor rather than contributing to the sample size, avoiding the inflation of degrees of freedom that arises when observational units are analyzed as independent participants."

**CLAIM-LOCK CHECK:** Supported. This is a conservative statistical-framing claim that constrains rather than expands inference, and it is consistent with the donor-aware modeling already used in the integrated analysis.

---

**SECTION:** Methods 2.1 / Limitations — Adolescent applicability of adult-derived preanalytical guidance

**PRIOR CLAIM:** Ex vivo fecal fermentation using stool from children and adolescents, including adolescents with obesity, is an established design; however, the stool storage, freeze–thaw, and inoculum-preservation literature is adult-derived.

**CITATIONS (max 3 lead + 2 backup):** Lead — Holmes 2020 (PMID 32788375); Mo 2025 (PMID 40895066). Backup — Aguirre 2015 (PMID 26222994); Gorzelak 2015 (PMID 26252519).

**CONFLICTS / NUANCE:** Holmes 2020 is the closest available precedent (adolescents with obesity, carbohydrate substrates, SCFA endpoints, 16S, donor-dependent responses), but it studied over-the-counter prebiotic supplements, not slow- vs rapid-digesting carbohydrate formulations. Preanalytical mechanisms are not known to be age-specific but have not been separately validated in adolescents.

**MANUSCRIPT-SAFE SENTENCE:** "Ex vivo fermentation of stool from children and adolescents, including adolescents with obesity, has been reported previously; specimen-handling procedures followed guidance derived from adult cohorts, as age-specific preanalytical validation for pediatric stool is not available."

**CLAIM-LOCK CHECK:** Supported as written. Do not extend Holmes 2020 to support obesity-group equivalence, preserved fermentation capacity, or substrate ranking — that paper reports donor-dependent variation and an explicit *lack* of positive correlation between fecal SCFA concentration, microbiota SCFA production capacity, and obesity markers.

---

## Verification notes

- 21 of 22 citations carry a PubMed-retrieved PMID **and** DOI. Hurlbert 1984 (#11) is not indexed in PubMed (1984 ecology journal); its DOI `10.2307/1942661` was verified against the ESA/Wiley journal landing page (`esajournals.onlinelibrary.wiley.com/doi/10.2307/1942661`) and multiple institutional copies. It is optional — Lazic 2018 (#10) alone carries the pseudoreplication argument if a PubMed-indexed-only reference list is preferred.
- No item in this report is marked `⚠ UNVERIFIED`.
- Searches executed: stool storage/freeze–thaw/16S; preanalytical SOPs; in vitro batch fermentation protocols and standardization; fresh-vs-frozen and cryopreserved fecal inocula; SCFA analytical stability; STORMS/MIxS/MBQC reporting standards; pseudoreplication and unit-of-analysis; technical replicates and replicate-sampling variance decomposition; antibiotic exposure and recovery; medication and diet covariates; adolescent/pediatric ex vivo fermentation; colonic pH control. Sources: PubMed (bio-research MCP), Consensus, and one targeted web verification.
