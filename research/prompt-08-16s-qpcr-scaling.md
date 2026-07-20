# Prompt 08 — 16S V3–V4 / DADA2 methods citations, genus-level functional limits, and qPCR-scaled abundance inference

**Manuscript target:** Methods §2.5–2.6; Discussion §4.6, §4.8 · **Citation job:** supply canonical methods citations for the amplicon + qPCR-scaling pipeline, and defensible caveat citations for what genus-level 16S and qPCR-scaled abundances can and cannot support · **Date:** 2026-07-19

All PMIDs/DOIs below were retrieved from PubMed (via the bio-research PubMed MCP tool) or Consensus during this run. Nothing is reconstructed from memory.

---

## Prioritized citation table

### A. Amplicon pipeline (Methods §2.5)

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 1 | DADA2: High-resolution sample inference from Illumina amplicon data | 2016 | Nat Methods | 10.1038/nmeth.3869 · PMID 27214047 | Canonical citation for ASV inference, error modeling, denoising. Mandatory in §2.5. |
| 2 | Exact sequence variants should replace operational taxonomic units in marker-gene data analysis | 2017 | ISME J | 10.1038/ismej.2017.119 · PMID 28731476 | Justifies ASVs (not OTUs) as the reporting unit; supports reusability/reproducibility rationale. |
| 3 | Evaluation of general 16S rRNA gene PCR primers for classical and next-generation sequencing-based diversity studies | 2013 | Nucleic Acids Res | 10.1093/nar/gks808 · PMID 22933715 | Source of the 341F/785R V3–V4 primer pair (S-D-Bact-0341-b-S-17 / S-D-Bact-0785-a-A-21); cite for primer choice and known coverage limits. |
| 4 | The SILVA ribosomal RNA gene database project: improved data processing and web-based tools | 2013 | Nucleic Acids Res | 10.1093/nar/gks1219 · PMID 23193283 | Reference database for taxonomy assignment (cite the release actually used). |
| 5 | Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy | 2007 | Appl Environ Microbiol | 10.1128/AEM.00062-07 · PMID 17586664 | The RDP naive-Bayes algorithm underlying `assignTaxonomy()` in DADA2; also documents that genus-level accuracy degrades on short variable-region fragments. |
| 6 | phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data | 2013 | PLoS ONE | 10.1371/journal.pone.0061217 · PMID 23630581 | Object handling, filtering, genus agglomeration (`tax_glom`), ordination. |
| 7 | Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data | 2018 | Microbiome | 10.1186/s40168-018-0605-2 · PMID 30558668 | decontam — the standard citation for negative-control / frequency-based contaminant filtering. |
| 8 | Towards standards for human fecal sample processing in metagenomic studies | 2017 | Nat Biotechnol | 10.1038/nbt.3960 · PMID 28967887 | DNA extraction is the single largest source of technical variation in fecal profiling; supports reporting extraction protocol explicitly. |
| 9 | Consistent and correctable bias in metagenomic sequencing experiments | 2019 | eLife | 10.7554/eLife.46923 · PMID 31502536 | Formal model of multiplicative taxon-specific measurement bias; grounds "within-study comparisons are valid, cross-protocol absolute values are not." |
| 10 | GTDB-Tk v2: memory friendly classification with the genome taxonomy database | 2022 | Bioinformatics | 10.1093/bioinformatics/btac672 · PMID 36218463 | **Use only if GTDB taxonomy was actually applied.** GTDB-Tk is genome-based; if the pipeline used SILVA only, cite #4/#5 and omit this. |

### B. Sequencing-depth handling (Methods §2.5)

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 11 | Waste not, want not: why rarefying microbiome data is inadmissible | 2014 | PLoS Comput Biol | 10.1371/journal.pcbi.1003531 · PMID 24699258 | The argument *against* rarefying for differential abundance; cite as one side of the depth controversy. |
| 12 | Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses | 2024 | mSphere | 10.1128/msphere.00354-23 · PMID 38251877 | The counter-position: rarefaction is preferred for α/β diversity. Pair with #11 to justify "rarefaction for diversity, un-rarefied counts for DA modeling." |
| 13 | Waste not, want not: revisiting the analysis that called into question the practice of rarefaction | 2023 | mSphere | 10.1128/msphere.00355-23 · PMID 38054712 | Companion reanalysis; useful if a reviewer challenges the depth-handling choice. |

### C. Compositionality and differential-abundance estimands (Methods §2.6; Discussion §4.6)

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 14 | Microbiome Datasets Are Compositional: And This Is Not Optional | 2017 | Front Microbiol | 10.3389/fmicb.2017.02224 · PMID 29187837 | The standard framing citation for why relative abundances carry no scale information. |
| 15 | Analysis of compositions of microbiomes with bias correction (ANCOM-BC) | 2020 | Nat Commun | 10.1038/s41467-020-17041-7 · PMID 32665548 | Sampling-fraction bias correction; the parent method. |
| 16 | Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures (ANCOM-BC2) | 2024 | Nat Methods | 10.1038/s41592-023-02092-7 · PMID 38158428 | **Primary ANCOM-BC2 citation.** Directly supports donor-as-repeated-measure designs — matches our paired donor × carbohydrate structure. |
| 17 | MaAsLin 3: refining and extending generalized multivariable linear models for meta-omic association discovery | 2026 | Nat Methods | 10.1038/s41592-025-02923-9 · PMID 41540124 | **Primary MaAsLin 3 citation.** Explicitly supports accounting for compositionality *experimentally* via qPCR or spike-ins, and separates abundance from prevalence associations. This is the citation that licenses the qPCR-scaled analysis. |
| 18 | Scale reliant mixed effects models enhance microbiome data analysis | 2026 | Microbiome | 10.1186/s40168-026-02377-x · PMID 41882807 | **Primary ALDEx3 citation.** Describes SR-MEM with user-specified scale uncertainty and random effects; the abstract names the ALDEx3 R package. |
| 19 | Addressing erroneous scale assumptions in microbe and gene set enrichment analysis | 2023 | PLoS Comput Biol | 10.1371/journal.pcbi.1011659 · PMID 37983251 | Foundational scale-reliant-inference paper: small errors in implicit scale assumptions can drive PPV as low as 9%. Best single citation for "why we modeled scale uncertainty rather than assuming a normalization." |
| 20 | Establishing microbial composition measurement standards with reference frames | 2019 | Nat Commun | 10.1038/s41467-019-10656-5 · PMID 31222023 | Reference frames; shows how naive cross-sample relative-abundance comparison generates false positives. |
| 21 | ALDEx2 (Fernandes et al., unified compositional tool) | 2014 | Microbiome | PMID 24910773 · ⚠ DOI not captured in this run — verify before use | Ancestor method of ALDEx3; cite only if ALDEx2 output is reported. Retrieve DOI at reference-management time. |
| 22 | Microbiome differential abundance methods produce different results across 38 datasets | 2022 | Nat Commun | 10.1038/s41467-022-28034-z · PMID 35039521 | Empirical basis for using more than one DA method — but see Conflicts below: the authors recommend a *consensus* approach, which is not our framing. |

### D. qPCR scaling / absolute abundance (Methods §2.6; Discussion §4.6, §4.8)

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 23 | Quantitative PCR provides a simple and accessible method for quantitative microbiota profiling | 2020 | PLoS ONE | 10.1371/journal.pone.0227285 · PMID 31940382 | **The direct methods precedent for our design:** total-bacterial qPCR run in parallel with library prep, used to scale NGS relative abundances into estimated absolute abundances. Cite in §2.6 for the scaling procedure itself. |
| 24 | Quantitative microbiome profiling links gut community variation to microbial load | 2017 | Nature | 10.1038/nature24460 · PMID 29143816 | Landmark QMP paper; up to 10-fold microbial-load variation among healthy individuals, and demonstrates that a relative-abundance "trade-off" can be a pure compositional artifact. |
| 25 | Absolute quantification of microbial taxon abundances | 2017 | ISME J | 10.1038/ismej.2016.117 · PMID 27612291 | **The cleanest citation for the enrichment-vs-outgrowth distinction:** "enrichment of taxa (increase in relative abundance) does not necessarily relate to the outgrowth of taxa (increase in absolute abundance)," and concludes both estimands should be reported. |
| 26 | Relative and Quantitative Rhizosphere Microbiome Profiling Results in Distinct Abundance Patterns | 2022 | Front Microbiol | 10.3389/fmicb.2021.798023 · PMID 35140695 | Explicit precedent for reporting relative and qPCR-estimated absolute abundance side by side as *complementary* information, including cases where the two disagree in direction. |
| 27 | Fecal microbial load is a major determinant of gut microbiome variation and a confounder for disease associations | 2024 | Cell | 10.1016/j.cell.2024.10.022 · PMID 39541968 | Microbial load, not the condition, explains much apparent disease-associated compositional change. Strong Discussion §4.8 caveat citation. |
| 28 | Uncertainty Modeling Outperforms Machine Learning for Microbiome Data Analysis | 2025 (preprint) | bioRxiv | 10.1101/2025.09.12.675956 · PMID 41000811 | ⚠ **Preprint — not peer-reviewed.** Counterweight to #27: ML-imputed microbial load fails to generalize; propagating scale uncertainty is preferable. Use only as a supporting citation. |

### E. Genus-level 16S and functional attribution (Discussion §4.6)

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 29 | Correcting for 16S rRNA gene copy numbers in microbiome surveys remains an unsolved problem | 2018 | Microbiome | 10.1186/s40168-018-0420-9 · PMID 29482646 | **The key terminology citation.** GCN prediction tools explain <10% of variance in some cases; authors recommend *against* routine copy-number correction. This is why our scaled quantity must be called 16S *gene-copy* equivalents, not cell counts. |
| 30 | rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea | 2015 | Nucleic Acids Res | 10.1093/nar/gku1201 · PMID 25414355 | Documents that rrn copy number varies across the tree and "introduces systematic bias when measuring community composition." |
| 31 | Evaluation of 16S rRNA gene sequencing for species and strain-level microbiome analysis | 2019 | Nat Commun | 10.1038/s41467-019-13036-1 · PMID 31695033 | Short-read variable-region targeting *cannot* achieve species/strain resolution; also documents intragenomic 16S copy variation. Directly limits how far genus-level calls can be pushed. |
| 32 | Function and functional redundancy in microbial systems | 2018 | Nat Ecol Evol | 10.1038/s41559-018-0519-1 · PMID 29662222 | Taxonomically distinct organisms encode the same energy-yielding functions, and the taxa carrying a function vary with little effect on the function. Core citation for "taxon shift ≠ function shift." |
| 33 | Understanding the effects of diet on bacterial metabolism in the large intestine | 2007 | J Appl Microbiol | 10.1111/j.1365-2672.2007.03322.x · PMID 17448155 | Butyrate producers are distributed across clostridial clusters IV and XIVa, and metabolic output depends on substrate *and* environment (pH). Supports "a genus label does not determine fermentation output." |
| 34 | PICRUSt2 for prediction of metagenome functions | 2020 | Nat Biotechnol | 10.1038/s41587-020-0548-6 · PMID 32483366 | Cite **only** if predicted-function output is reported; otherwise cite in the Discussion as the approach we deliberately did *not* use to attribute metabolic function. |

---

## Prior-claim bullets (mapped to citation numbers)

- **ASV inference from V3–V4 reads with DADA2 is standard practice**, and ASVs are preferred over OTUs as the analysis unit because they are database-independent, reproducible labels. [1, 2]
- **The 341F/785R V3–V4 pair is an evaluated, coverage-benchmarked choice**, selected in silico against SILVA — but every primer pair carries phylum-level coverage bias. [3, 4]
- **Naive-Bayes taxonomy assignment against SILVA is the default DADA2 path**, and its own validation shows genus-level error rates rise for short variable-region fragments. [4, 5]
- **Contaminant removal via negative controls (decontam) and prevalence/depth filtering are expected reporting elements**, not optional extras. [6, 7]
- **Fecal DNA extraction protocol dominates technical variation**, and taxon-specific measurement bias is multiplicative and protocol-specific — so within-study contrasts are interpretable while absolute cross-study values are not. [8, 9]
- **Depth handling is genuinely contested.** Rarefying was argued to be inadmissible for differential abundance [11]; more recent simulation work argues rarefaction is the most robust control for uneven effort in α/β diversity [12, 13]. A split strategy (rarefaction for diversity, un-rarefied counts with model-based scaling for DA) is defensible if stated.
- **Relative abundances are compositional and carry no scale information**; conclusions drawn from them are conditional on an implicit assumption about total load. [14, 19, 20]
- **ANCOM-BC2 is the appropriate relative-abundance estimator for a repeated-measures donor design**, having been developed for covariate adjustment and repeated measures. [15, 16]
- **MaAsLin 3 explicitly supports experimentally-informed compositionality correction via qPCR or spike-ins**, which is precisely the qPCR-scaled estimand we report. [17]
- **ALDEx3 / SR-MEM propagates user-specified uncertainty about the unmeasured scale through a mixed-effects model** rather than fixing scale by normalization, and can incorporate external qPCR/flow measurements. [18, 19]
- **Total-bacterial 16S qPCR run in parallel with library prep is an established, accessible route to estimated absolute abundances.** [23, 24]
- **A rise in relative abundance is not evidence of outgrowth.** Props et al. state this directly and conclude that both relative and absolute abundances should be reported for a comprehensive interpretation. [25]
- **Relative and qPCR-estimated absolute profiles can disagree in direction**, and reporting both is an established practice that adds information rather than adjudicating a winner. [24, 25, 26]
- **Microbial load is itself a major driver of apparent compositional differences** and can confound group contrasts if unmeasured. [24, 27]
- **16S gene copy number varies across taxa and cannot be reliably predicted**, so a qPCR-scaled quantity is an estimate of 16S gene-copy (target) equivalents, not of cells. [29, 30, 31]
- **Short-read V3–V4 sequencing does not resolve species or strains**, and intragenomic 16S copies differ within a single genome. [31]
- **Metabolic function is functionally redundant across taxa and environmentally contingent**, so genus-level abundance change does not license a mechanistic metabolic attribution. [32, 33]

---

## Conflicts / nuance

1. **Rarefaction: an unresolved methodological dispute.** McMurdie & Holmes [11] and Schloss [12, 13] reach opposite conclusions. Do not present either as settled. State the choice, cite both sides, and note that the DA models used here (ANCOM-BC2, MaAsLin 3, ALDEx3) operate on un-rarefied counts by design.

2. **"Consensus approach" vs. "complementary estimands" — an important framing conflict.** Nearing et al. [22] explicitly recommend "a consensus approach based on multiple differential abundance methods." That is method-voting, and it is *not* our stated rationale. Our three tools do not estimate the same quantity: ANCOM-BC2 targets a relative/compositional contrast, MaAsLin 3 (qPCR-informed) targets a scaled-abundance contrast, and ALDEx3 targets a contrast that is deliberately partially identified under scale uncertainty. Cite [22] only for the empirical fact that DA methods disagree; cite [19, 25, 26] for the framing that the disagreement is partly *definitional* rather than a reliability contest. Agreement across the three is reassuring but is not proof of truth, and disagreement is not automatically an error — it may reflect a real difference between relative and absolute change.

3. **qPCR scaling reduces one bias but introduces others.** Total-bacterial qPCR is itself subject to universal-primer coverage bias, copy-number variation, extraction efficiency differences, and inhibitors. Jian et al. [23] present the method as accessible and useful while noting challenges; Louca et al. [29] and Stoddard et al. [30] establish that the copy-number problem is unresolved. Net effect: the scaled quantity is an *estimate* with its own error, not a ground truth against which the relative analysis is scored.

4. **Load prediction from sequence data is contested.** Nishijima et al. [27] used a machine-learning load predictor at scale; Konnaris et al. [28, preprint] report that such models fail to generalize across studies. We measured load rather than predicting it, so this dispute does not bear on our estimates directly — but it argues against ever treating a scaled abundance as though the scale were known exactly.

5. **PICRUSt2 [34] is not a substitute for measured function.** If it is not run, say so; if it is run, its output is a prediction constrained by nearest-sequenced-taxon distance and should not be described as observed metabolic capacity. Louca et al. [29] specifically evaluated PICRUSt-family copy-number prediction accuracy and found it poor across many taxa.

6. **Genus agglomeration loses information asymmetrically.** Aggregating ASVs to genus can merge organisms with different substrate preferences [31, 32, 33]. This weakens, not strengthens, any inference from a genus-level change to a fermentation outcome.

---

## What we must not overclaim (given claim lock)

- **Do not call qPCR-scaled values "cell counts," "absolute abundance," or "bacterial numbers" without qualification.** Copy-number variation is unpredictable [29, 30] and 16S copies vary even within one genome [31]. Approved terms: *qPCR-scaled genus abundance*, *estimated 16S rRNA gene-copy equivalents*, *estimated 16S target equivalents*. Every figure axis and table header should carry the qualifier at least once.
- **Do not report absolute expansion of any genus (including *Fusicatenibacter*) as demonstrated.** A relative-abundance increase is compatible with no change or even a decline in scaled abundance [24, 25]. If the scaled analysis does not show a directionally concordant, well-supported increase, the correct statement is that relative enrichment was observed without evidence of absolute expansion.
- **Do not attribute SCFA production to any genus from 16S data.** Butyrate production is polyphyletic and environment-dependent [32, 33], and V3–V4 short reads cannot resolve the species or strains that would carry the relevant pathways [31]. No citation in this set supports a genus-to-metabolite mechanistic link in this design.
- **Do not frame agreement among ANCOM-BC2 / MaAsLin 3 / ALDEx3 as validation.** They estimate different quantities under different scale assumptions [17, 18, 19]. "All three agreed" is a robustness observation, not evidence of a true effect.
- **Do not describe qPCR scaling as removing compositionality.** It supplies an estimate of scale with its own uncertainty [19, 23, 28]; it does not convert the data into an error-free absolute measurement.
- **Do not use PICRUSt2 or any predicted-function output to claim metabolic capacity** [29, 34].
- **Do not extend any of these citations to obesity-group equivalence, responder phenotypes, SDC-over-RDC superiority, or dietary-response prediction.** Nothing retrieved here speaks to those claims, and none should be recruited to support them.

---

## Methods citation checklist

Copy-paste checklist for Methods §2.5–2.6. Mark each line as **cited / N-A / verify-version**.

**§2.5 — Amplicon sequencing and bioinformatics**

- [ ] DNA extraction protocol named; technical-variation caveat → Costea 2017 [8]
- [ ] V3–V4 primers identified by name/sequence → Klindworth 2013 [3]
- [ ] Read processing, filtering, denoising, merging, chimera removal → Callahan 2016 (DADA2) [1]
- [ ] ASVs stated as the analysis unit (not OTUs) → Callahan 2017 [2]
- [ ] Taxonomy assignment method → Wang 2007 (naive Bayes) [5]
- [ ] Reference database **with release/version number** → Quast 2013 (SILVA) [4]; GTDB-Tk [10] only if GTDB used
- [ ] Negative-control-based contaminant removal → Davis 2018 (decontam) [7]
- [ ] Prevalence / minimum-count / depth filters — state exact thresholds; no separate citation required
- [ ] Genus agglomeration and data handling → McMurdie & Holmes 2013 (phyloseq) [6]
- [ ] Depth handling stated explicitly (rarefied for diversity? un-rarefied for DA?) → [11] + [12] (and [13] if challenged)
- [ ] Measurement-bias framing for within- vs cross-study comparison → McLaren 2019 [9]

**§2.6 — Total bacterial qPCR and differential abundance**

- [ ] Total bacterial 16S qPCR assay, standard curve, and units reported
- [ ] Parallel qPCR + sequencing scaling procedure → Jian 2020 [23]
- [ ] Scaled quantity named as *estimated 16S gene-copy equivalents* (never "cell counts") → Louca 2018 [29], Stoddard 2015 [30]
- [ ] Compositionality rationale → Gloor 2017 [14]
- [ ] Relative-abundance DA model, donor as repeated measure → Lin & Peddada 2024 (ANCOM-BC2) [16] (+ [15] for the parent method)
- [ ] qPCR-informed DA model → Nickols 2026 (MaAsLin 3) [17]
- [ ] Scale-uncertainty DA model → McGovern & Silverman 2026 (ALDEx3 / SR-MEM) [18], with the scale-assumption rationale from McGovern 2023 [19]
- [ ] ALDEx2 [21] cited **only** if ALDEx2 results are shown — ⚠ retrieve DOI before use
- [ ] Multiplicity control (FDR method and level) stated per model
- [ ] Explicit statement that the three models estimate **different estimands** and are reported side by side, not vote-counted → [19, 25, 26]
- [ ] Statement that all 16S/qPCR analyses are secondary/exploratory relative to the SCFA primary endpoints

---

## Downstream synthesis blocks

### Block 1

**SECTION:** Methods §2.5 — amplicon processing
**PRIOR CLAIM:** ASV-level inference from V3–V4 amplicons with DADA2, SILVA taxonomy assignment, negative-control contaminant filtering, and genus agglomeration constitute standard, reproducible practice for community profiling.
**CITATIONS (lead):** Callahan 2016 [1]; Callahan 2017 [2]; Davis 2018 (decontam) [7] — **backup:** Quast 2013 (SILVA) [4]; McMurdie & Holmes 2013 (phyloseq) [6]
**CONFLICTS / NUANCE:** Genus-level assignment accuracy from short variable-region fragments is imperfect [5]; primer coverage is not uniform across phyla [3]; extraction protocol drives substantial technical variation [8].
**MANUSCRIPT-SAFE SENTENCE:** "Paired-end V3–V4 reads were denoised into amplicon sequence variants with DADA2, taxonomy was assigned by naive Bayesian classification against SILVA (release X), putative contaminants were removed using negative controls in decontam, and ASVs were agglomerated to genus in phyloseq."
**CLAIM-LOCK CHECK:** Supported — descriptive methods reporting only.

### Block 2

**SECTION:** Methods §2.6 — qPCR scaling
**PRIOR CLAIM:** Running total-bacterial 16S qPCR in parallel with library preparation and multiplying genus relative abundances by total load yields an estimate of absolute abundance that mitigates compositional constraints.
**CITATIONS (lead):** Jian 2020 [23]; Vandeputte 2017 [24]; Gloor 2017 [14] — **backup:** Props 2017 [25]; McGovern 2023 [19]
**CONFLICTS / NUANCE:** The scaled quantity inherits qPCR primer-coverage bias, extraction-efficiency differences, and — decisively — unresolved 16S copy-number variation [29, 30, 31]. It is an estimate carrying its own uncertainty, not a ground-truth measurement.
**MANUSCRIPT-SAFE SENTENCE:** "Genus relative abundances were multiplied by total bacterial 16S rRNA gene copies quantified by qPCR on the same DNA extracts, yielding qPCR-scaled genus abundances expressed as estimated 16S rRNA gene-copy equivalents; because 16S copy number varies among taxa and cannot be reliably predicted, these values are not interpreted as cell counts."
**CLAIM-LOCK CHECK:** Supported, provided the "gene-copy equivalents" qualifier is retained wherever the quantity is named.

### Block 3

**SECTION:** Methods §2.6 — differential abundance strategy
**PRIOR CLAIM:** ANCOM-BC2, MaAsLin 3, and ALDEx3 were applied because they estimate different, complementary quantities under different assumptions about the unmeasured biological scale — not to adjudicate a single answer by majority.
**CITATIONS (lead):** Lin & Peddada 2024 (ANCOM-BC2) [16]; Nickols 2026 (MaAsLin 3) [17]; McGovern & Silverman 2026 (ALDEx3) [18] — **backup:** McGovern 2023 [19]; Nearing 2022 [22]
**CONFLICTS / NUANCE:** Nearing et al. [22] recommend a *consensus* approach across DA methods; we deliberately do not adopt method-voting, because the three estimands are not interchangeable. Cite [22] for the empirical fact of method disagreement only.
**MANUSCRIPT-SAFE SENTENCE:** "Differential abundance was assessed with three complementary models rather than a single method: ANCOM-BC2 for the relative-abundance contrast with donor as a repeated measure, MaAsLin 3 using the qPCR-derived total load to account for compositionality experimentally, and ALDEx3 to propagate uncertainty in the unmeasured scale; because these approaches target different estimands, results are reported alongside one another rather than combined into a consensus call."
**CLAIM-LOCK CHECK:** Supported — but the sentence must not be paraphrased into "concordance across methods confirms the finding."

### Block 4

**SECTION:** Discussion §4.6 — limits of genus-level 16S for functional attribution
**PRIOR CLAIM:** Genus-level 16S profiling cannot establish which organisms generated the observed SCFA production, because metabolic functions are redundant across taxa, environmentally contingent, and encoded below the resolution of short-read V3–V4 sequencing.
**CITATIONS (lead):** Louca 2018 (functional redundancy) [32]; Johnson 2019 [31]; Louis 2007 [33] — **backup:** Louca 2018 (copy number) [29]; Douglas 2020 (PICRUSt2) [34]
**CONFLICTS / NUANCE:** Predicted-function tools exist [34] but their underlying copy-number/gene predictions are unreliable for taxa distant from sequenced genomes [29]; they would not resolve this limitation.
**MANUSCRIPT-SAFE SENTENCE:** "Because short-read V3–V4 sequencing does not resolve species or strains, and because fermentative capacity is distributed redundantly across taxa and modulated by substrate and environmental conditions, the genus-level compositional shifts observed here cannot be used to attribute the measured SCFA production to particular organisms."
**CLAIM-LOCK CHECK:** Supported as a limitation. Any inversion of this sentence into a mechanistic attribution (e.g. naming a genus as the butyrate source) would be **unsupported**.

### Block 5

**SECTION:** Discussion §4.8 — relative increase without absolute expansion
**PRIOR CLAIM:** An increase in a genus's relative abundance can occur without any increase in its absolute abundance, because relative abundances are constrained to sum to one and total microbial load varies substantially.
**CITATIONS (lead):** Props 2017 [25]; Vandeputte 2017 [24]; Gloor 2017 [14] — **backup:** Morton 2019 [20]; Nishijima 2024 [27]
**CONFLICTS / NUANCE:** The converse also holds — the qPCR-scaled estimate is itself uncertain [19, 29], so a null scaled result is not positive proof of no expansion; it is absence of evidence for expansion. State it that way.
**MANUSCRIPT-SAFE SENTENCE:** "Because sequencing yields proportions rather than counts, enrichment in relative abundance need not correspond to outgrowth in absolute terms; in this dataset the qPCR-scaled estimates did not provide evidence of absolute expansion, and given the uncertainty inherent in gene-copy-based scaling we treat this as an absence of supporting evidence rather than as demonstrated stability."
**CLAIM-LOCK CHECK:** Supported. Guards directly against the locked-out claim of absolute *Fusicatenibacter* expansion, while also avoiding the mirror-image overclaim of proven non-expansion.

### Block 6

**SECTION:** Discussion §4.8 — microbial load as a variable in its own right
**PRIOR CLAIM:** Total microbial load varies substantially between individuals and can itself drive apparent compositional differences, so load should be measured and reported rather than assumed constant.
**CITATIONS (lead):** Vandeputte 2017 [24]; Nishijima 2024 [27] — **backup:** Morton 2019 [20]; Konnaris 2025 preprint [28]
**CONFLICTS / NUANCE:** Nishijima et al. inferred load with a machine-learning predictor; Konnaris et al. [28] (preprint, not peer-reviewed) report that such predictors generalize poorly. We measured load by qPCR, which avoids that dispute but does not escape assay-level bias [23, 29].
**MANUSCRIPT-SAFE SENTENCE:** "Total bacterial load varies markedly among donors and can itself account for apparent differences in community composition; we therefore measured load directly by qPCR in every sample rather than assuming it constant across conditions."
**CLAIM-LOCK CHECK:** Supported — descriptive and design-justifying; makes no group-level or clinical claim.

---

## Verification notes

- Items [1]–[20], [22]–[27], [29]–[34] were retrieved with full metadata (title, journal, year, DOI, PMID) from PubMed during this run.
- **[21] ALDEx2 (Fernandes et al. 2014, Microbiome, PMID 24910773)** — PMID confirmed by PubMed search; full metadata/DOI was not fetched. `⚠ UNVERIFIED DOI — needs manual check` before it enters the reference manager. Cite only if ALDEx2 results are actually reported.
- **[28] Konnaris et al. 2025 (bioRxiv, PMID 41000811, DOI 10.1101/2025.09.12.675956)** — verified, but is a **preprint and not peer-reviewed**. Flag as such if cited, or drop.
- **[10] GTDB-Tk v2** — verified, but is a genome-based classifier. Include only if GTDB taxonomy was genuinely used in this pipeline; otherwise omit to avoid a misleading methods citation.
- **[16] vs [37205444]** — the ANCOM-BC2 Research Square preprint (PMID 37205444) exists; always cite the peer-reviewed Nat Methods version [16] instead.
- **[17] vs [39713460]** — likewise, the MaAsLin 3 bioRxiv preprint exists; cite the peer-reviewed Nat Methods version [17].
- SILVA and any database release numbers must be filled in from the actual pipeline configuration before submission.
