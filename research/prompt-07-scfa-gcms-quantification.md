# Prompt 07 — SCFA quantification by PFBBr/GC-MS (NCI): methods citations, QC, and reporting language

**Manuscript target:** Methods §2.4 (SCFA analytical chemistry); Results/Discussion language for SCFA endpoints · **Citation job:** Supply facility-confirmable, peer-reviewed analytical-chemistry citations for a PFBBr-derivatization GC-MS/NCI SCFA workflow (plus acceptable alternatives), reporting-standard citations for units/dilution/below-LOQ handling, citations motivating succinate and 5-aminovalerate as related fermentation metabolites, and the citation basis for using *concentration / net accumulation / net output* rather than *production* · **Date:** 2026-07-19

> **Verification note.** Every PMID and DOI below was retrieved from PubMed via tool query (per `_SHARED_CONTEXT.md` citation-integrity rules). Nothing was reconstructed from memory. One item that surfaced only as a conference abstract via Consensus, and could not be confirmed in PubMed, is flagged `⚠ UNVERIFIED`. No vendor application notes or product literature are included, per the prompt constraint.

---

## Prioritized citation table

### A. PFBBr / negative-chemical-ionization core methodology

| # | Title | Year | Journal | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 1 | Isotopomer enrichment assay for very short chain fatty acids and its metabolic applications | 2011 | Anal Biochem | PMID 21112315 · [10.1016/j.ab.2010.11.030](https://doi.org/10.1016/j.ab.2010.11.030) | **Lead canonical citation.** The closest published analogue to our workflow: PFBBr derivatization of formic, acetic, propionic, butyric and pentanoic acids, detected by GC-MS in negative chemical ionization. Explicitly optimizes pH, temperature, reaction time and PFBBr:sample ratio, and reports precision, stability and accuracy validation. Cite for the derivatization chemistry and the NCI rationale. |
| 2 | Simultaneous determination of trace levels of nine haloacetic acids in biological samples as their pentafluorobenzyl derivatives by GC/MS-MS in electron capture negative ion chemical ionization mode | 2003 | Anal Chem | PMID 14632119 · [10.1021/ac034036w](https://doi.org/10.1021/ac034036w) | Foundational demonstration that PFBBr esters of small carboxylic acids ionize efficiently under ECNICI to give the carboxylate anion, with pg/mL LODs in plasma and urine. Cite for why ECNICI is chosen and for the LOD/accuracy/precision framework (intra- vs inter-day accuracy and CV reported separately). |
| 3 | GC-MS and GC-MS/MS measurement of ibuprofen … after ethyl acetate extraction and pentafluorobenzyl bromide derivatization | 2016 | J Chromatogr B | PMID 27343144 · [10.1016/j.jchromb.2016.06.014](https://doi.org/10.1016/j.jchromb.2016.06.014) | Modern, fully validated PFBBr/ECNICI carboxylic-acid assay using a stable-isotope-labelled internal standard and acetonitrile/DIPEA derivatization. Cite for the internal-standard and base-catalyst conditions and for the `[M−PFB]⁻` carboxylate-anion monitoring scheme. |
| 4 | Determination of 2-, 3-, 4-methylpentanoic and cyclohexanecarboxylic acids … SPE and gas chromatography-negative chemical ionization mass spectrometry | 2015 | J Chromatogr A | PMID 25601317 · [10.1016/j.chroma.2014.12.074](https://doi.org/10.1016/j.chroma.2014.12.074) | PFBBr room-temperature derivatization (30 min) of short/branched carboxylic acids with GC-NCI-MS; reports detection limits, linear range, reproducibility (RSD < 10%) and ~100% signal recovery. Useful supporting citation for mild derivatization conditions and for branched-chain analytes. |

### B. Acceptable alternative GC-MS SCFA methods commonly used in microbiome studies

| # | Title | Year | Journal | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 5 | A targeted metabolomic protocol for short-chain fatty acids and branched-chain amino acids | 2013 | Metabolomics | PMID 23997757 · [10.1007/s11306-013-0500-6](https://doi.org/10.1007/s11306-013-0500-6) | The propyl chloroformate (PCF) aqueous-derivatization alternative, widely adopted in microbiome/obesity work. Reports intra- and inter-day RSD < 10% and stability at room temperature and −20 °C. Cite when justifying derivatization choice against alternatives. |
| 6 | A sensitive GC/MS detection method for analyzing microbial metabolites short chain fatty acids in fecal and serum samples | 2019 | Talanta | PMID 30683360 · [10.1016/j.talanta.2018.12.049](https://doi.org/10.1016/j.talanta.2018.12.049) | The BSTFA/TMS silylation alternative with sodium-sulfate dehydration; LODs 0.064–0.067 µM, R² > 0.999, RSD < 2%. Heavily cited benchmark for what "good" SCFA GC-MS performance looks like in fecal matrices. |
| 7 | A Gas Chromatography Mass Spectrometry-Based Method for the Quantification of Short Chain Fatty Acids | 2022 | Metabolites | PMID 35208244 · [10.3390/metabo12020170](https://doi.org/10.3390/metabo12020170) | Derivatization-free GC-MS across plasma, feces, cecum, liver and adipose. **Most useful for our QC language**: quantifies using deuterated internal standards *plus* external calibration curves, and explicitly tests and rules out matrix effects per matrix. Cite for the dual IS + external-calibration design and for matrix-effect assessment. |
| 8 | Alternative method for gas chromatography-mass spectrometry analysis of short-chain fatty acids in faecal samples | 2012 | J Sep Sci | PMID 22865755 · [10.1002/jssc.201101121](https://doi.org/10.1002/jssc.201101121) | Direct-injection fecal SCFA GC-MS using 4-methylvaleric acid as internal standard; reports recovery 65–105%, no matrix effects, LODs 0.49–4.31 µM. Cite as the standard non-isotopic internal-standard precedent in fecal work. |
| 9 | Orthogonal assessment of six extraction methods for GC-MS based short-chain fatty acid quantification in fecal samples | 2025 | J Chromatogr A | PMID 40527045 · [10.1016/j.chroma.2025.466132](https://doi.org/10.1016/j.chroma.2025.466132) | Recent head-to-head comparison of six extraction chemistries showing **no single method is optimal for all SCFAs** — recovery and sensitivity are analyte-specific. Cite to justify a method choice honestly rather than claiming universal superiority, and to motivate per-analyte recovery reporting. |

### C. Interpreting culture supernatant vs. fecal SCFA measurements

| # | Title | Year | Journal | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 10 | Pitfalls in short-chain fatty acid research: A methodological review | 2019 | Anim Sci J | PMID 30370625 · [10.1111/asj.13118](https://doi.org/10.1111/asj.13118) | **Lead language-guardrail citation.** States directly that lumen or fecal SCFA concentration does *not* reflect rate of production, and that such parameters "should not be used as measures of SCFA production or absorption." Also warns against over-emphasis on n-butyrate. This is the single most quotable source for our concentration-vs-production wording. |
| 11 | Short-chain fatty acid kinetics and concentrations are higher after inulin supplementation in young and older adults: a randomized trial | 2025 | Am J Clin Nutr | PMID 40274191 · [10.1016/j.ajcnut.2025.04.018](https://doi.org/10.1016/j.ajcnut.2025.04.018) | Empirical human demonstration that **plasma but not fecal** SCFA concentrations correlate with tracer-measured production; fecal concentration changes and production changes dissociate. Strongest modern evidence that a concentration endpoint is not a production endpoint. |
| 12 | Systemic availability and metabolism of colonic-derived short-chain fatty acids in healthy subjects: a stable isotope study | 2017 | J Physiol | PMID 27510655 · [10.1113/JP272613](https://doi.org/10.1113/JP272613) | Quantifies the *in vivo* sinks absent from a closed ex vivo vessel: only 36%/9%/2% of colonic acetate/propionate/butyrate reach the systemic circulation, with 24% acetate→butyrate microbial interconversion. Cite to explain why an ex vivo supernatant measurement and an in vivo fecal measurement are not interchangeable, in either direction. |
| 13 | Key bacterial taxa and metabolic pathways affecting gut short-chain fatty acid profiles in early life | 2021 | ISME J | PMID 33723382 · [10.1038/s41396-021-00937-7](https://doi.org/10.1038/s41396-021-00937-7) | Longitudinal fecal organic-acid profiling in which succinate, lactate and formate are measured and interpreted **alongside** acetate/propionate/butyrate as phase-defining fermentation intermediates. Precedent for reporting a wider organic-acid panel rather than the three canonical SCFAs alone. |

### D. Reporting standards: units, dilution correction, below-LOQ handling, batch/QC

| # | Title | Year | Journal | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 14 | Mass spectrometry-based metabolomics: a guide for annotation, quantification and best reporting practices | 2021 | Nat Methods | PMID 34239102 · [10.1038/s41592-021-01197-1](https://doi.org/10.1038/s41592-021-01197-1) | **Lead reporting-standards citation.** Community guidelines covering sample preparation, replication and randomization, quantification, recovery and recombination, ion suppression and peak misidentification, for both LC-MS and GC-MS. Cite as the umbrella standard our Methods §2.4 conforms to. |
| 15 | Use cases, best practice and reporting standards for metabolomics in regulatory toxicology (MERIT) | 2019 | Nat Commun | PMID 31292445 · [10.1038/s41467-019-10900-y](https://doi.org/10.1038/s41467-019-10900-y) | Performance and reporting standards for **targeted metabolite data** specifically, from a multi-sector expert consortium. Cite for the requirement to state LOD/LOQ, calibration range, and QC acceptance criteria explicitly rather than by reference. |
| 16 | Development and Application of Ultra-Performance Liquid Chromatography-TOF MS for Precision Large Scale Urinary Metabolic Phenotyping | 2016 | Anal Chem | PMID 27479709 · [10.1021/acs.analchem.6b01481](https://doi.org/10.1021/acs.analchem.6b01481) | Canonical description of the **pooled-QC-interspersed** run design: repeated pooled QC injections distributed throughout each analytical batch, with retention-time and peak-area RSD as the acceptance metric, and QC-dilution-series filtering. Cite for our batch-handling and drift-monitoring design. |
| 17 | Effect of different pooled QC samples on data quality during an inter-batch experiment | 2025 | Anal Bioanal Chem | PMID 39557686 · [10.1007/s00216-024-05646-6](https://doi.org/10.1007/s00216-024-05646-6) | Shows that *how* the pooled QC is constituted (pooled pre- vs post-extraction) materially changes inter-batch correction and downstream feature selection. Cite to justify stating explicitly at which step our QC pool was formed. |

### E. Motivating succinate and 5-aminovalerate as related fermentation metabolites

| # | Title | Year | Journal | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 18 | Phylogenetic distribution of three pathways for propionate production within the human gut microbiota | 2014 | ISME J | PMID 24553467 · [10.1038/ismej.2014.14](https://doi.org/10.1038/ismej.2014.14) | **Lead succinate citation.** Establishes the succinate pathway as the dominant route to propionate from hexose sugars in the human colonic microbiota. This is the mechanistic reason succinate belongs in a carbohydrate-fermentation panel: it is the principal propionate precursor, and residual succinate is informative about incomplete conversion. |
| 19 | Fructo-oligosaccharides promote butyrate production over citrus pectin during in vitro fermentation by colonic inoculum | 2024 | Anaerobe | PMID 39393609 · [10.1016/j.anaerobe.2024.102919](https://doi.org/10.1016/j.anaerobe.2024.102919) | Direct **in vitro batch fermentation** precedent in which succinate is interpreted as a cross-feeding intermediate between succinate-producing and succinate-utilizing taxa, and substrate chemistry shifts the acetate/propionate/butyrate balance. Closest design analogue for reporting succinate in a carbohydrate-comparison fermentation. |
| 20 | Gut microbial activity, implications for health and disease: the potential role of metabolite analysis | 2012 | J Proteome Res | PMID 23116228 · [10.1021/pr300637d](https://doi.org/10.1021/pr300637d) | Review of **proteolytic** (amino-acid-derived) gut fermentation products and the analytical methods used to detect them. Cite as the general rationale for reporting amino-acid-derived metabolites such as 5-aminovalerate (a lysine-fermentation product) as a marker of the saccharolytic-versus-proteolytic balance when carbohydrate is limiting. |
| 21 | Microbe-Immune Crosstalk: Evidence That T Cells Influence the Development of the Brain Metabolome | 2022 | Int J Mol Sci | PMID 35328680 · [10.3390/ijms23063259](https://doi.org/10.3390/ijms23063259) | Concrete precedent measuring and reporting **5-aminovalerate** as a named gut/cecal microbial metabolite in a profiling panel (it is a listed keyword of the paper). Cite as evidence 5-aminovalerate is an established reported analyte, not an ad hoc addition. |
| 22 | Human Milk Oligosaccharide-Stimulated *Bifidobacterium* Species Contribute to Prevent Later Respiratory Tract Infections | 2021 | Microorganisms | PMID 34576834 · [10.3390/microorganisms9091939](https://doi.org/10.3390/microorganisms9091939) | Reports **succinate and 5-aminovalerate together with acetate, butyrate and propionate** in a single stool organic-acid panel in a human carbohydrate (HMO) intervention. The tightest single precedent for our exact analyte list. |
| 23 | A novel Mediterranean diet-inspired supplement ameliorates cognitive, microbial, and metabolic deficits in a mouse model of low-grade inflammation | 2024 | Gut Microbes | PMID 38835220 · [10.1080/19490976.2024.2363011](https://doi.org/10.1080/19490976.2024.2363011) | Additional recent precedent reporting 5-aminovalerate as a quantified microbiota-linked metabolite. Backup citation only. |

### Flagged item

| # | Title | Year | Journal | Status |
|---|---|---|---|---|
| 24 | Systemic but not fecal short-chain fatty acid (SCFA) concentrations reflect intestinal SCFA production as assessed by novel stable isotope method (Kirschner et al.) | 2023 | Clin Nutr ESPEN | **⚠ UNVERIFIED — needs manual check.** Surfaced via Consensus only; no PMID/DOI returned and not confirmable in PubMed. Appears to be a conference abstract. **Do not cite.** Use #11 (Kirschner 2025, Am J Clin Nutr, PMID 40274191) instead — it is the peer-reviewed full paper making the same point. |

---

## Methods citation block

*Drop-in prose for Methods §2.4. Bracketed numbers map to the table above; swap for the journal's citation style.*

> **2.4 Short-chain fatty acid quantification.** Fermentation supernatants were analyzed for acetate, propionate and butyrate, together with succinate and 5-aminovalerate, by gas chromatography–mass spectrometry following pentafluorobenzyl bromide (PFBBr) derivatization with detection under negative chemical ionization. Derivatization and NCI detection followed established procedures for very-short-chain and short-chain carboxylic acids [1, 2], in which PFBBr esters ionize under electron-capture negative-ion chemical ionization to yield the corresponding carboxylate anion `[M−PFB]⁻`, monitored by selected-ion monitoring [1–3]. Derivatization conditions (pH, temperature, reaction time, and PFBBr-to-sample ratio) were those optimized for this analyte class [1], with mild room-temperature derivatization as described for short and branched carboxylic acids [4].
>
> Quantification used stable-isotope-labelled internal standards added prior to extraction, combined with external calibration curves prepared in matrix-matched blank fermentation medium, following the dual internal-standard-plus-external-calibration design validated for SCFA quantification across biological matrices [7]; where an isotopic analogue was unavailable, a non-endogenous carboxylic acid internal standard was used as is standard in fecal SCFA work [8]. Matrix effects were assessed by comparing calibration slopes in neat solvent against matrix-matched standards [7], and analyte-specific recovery was determined by spiking at low, medium and high concentrations, since extraction efficiency for SCFAs is known to be analyte-dependent rather than uniform across a panel [9].
>
> Method performance was characterized by linearity across the calibration range, intra- and inter-day accuracy and precision, and limits of detection and quantification (LOD/LOQ), reported per analyte [2, 6, 8]. Samples were randomized across the run and analyzed in batches, with pooled quality-control samples — constituted by combining aliquots of study supernatants prior to extraction — injected at the start of and at regular intervals throughout each batch to condition the column, monitor drift and quantify measurement precision [16, 17]. Retention-time and peak-area relative standard deviations in the pooled QC injections were used as batch acceptance criteria [16]. Acquisition, quantification and QC reporting conform to community best-practice guidelines for mass-spectrometry-based metabolomics [14] and to reporting standards for targeted metabolite data [15].
>
> Concentrations are reported in µmol/L of fermentation supernatant, corrected for any dilution applied during sample preparation. Values below the LOQ were handled by a single pre-specified rule applied uniformly to all analytes and groups, and the number and distribution of such values are reported per analyte and per group [14, 15]. Because the incubation is a closed ex vivo system, endpoints are expressed as the **net change in concentration from 0 to 48 h** (net accumulation), which reflects the balance of microbial formation and consumption within the vessel rather than a production rate [10–12].

---

## Six reporting best-practice bullets, with citations

1. **Report units unambiguously, and state the denominator.** Give SCFA concentrations in molar units per volume of supernatant for culture work (µmol/L), and state explicitly whether any per-mass normalization was applied. Molar units make the analytes directly comparable across the panel and avoid the mass-basis ambiguity that plagues fecal reporting; per-analyte units and calibration ranges are an explicit expectation of targeted-metabolite reporting standards [15], within the broader quantification guidance of [14]. For our design, the supernatant volume is a defined, measured quantity — this is a genuine advantage over fecal work and should be stated as such.

2. **Report dilution correction as a stated arithmetic step, not an implicit one.** Any dilution applied to bring high-abundance acetate into the calibration range must be recorded and the correction factor reported, because dilution is applied unequally across analytes spanning a wide dynamic range. Guidance on quantification, recovery and recombination in complex mixtures treats dynamic-range handling as a primary source of quantitative error [14]; the same source flags ion suppression, which changes with dilution, as requiring explicit assessment.

3. **Pre-specify a single below-LOQ rule and report the censoring rate per analyte per group.** Choose one rule (for example, substitution at LOQ/2, or retention of the instrument-reported value with a censoring flag) before unblinding, apply it identically to every analyte and group, and report how many observations it affected in each group. Reporting standards for targeted metabolite data require LOD/LOQ to be stated and the treatment of sub-limit values to be transparent [15, 14]. This matters directly for us: 5-aminovalerate and succinate are plausibly near the LOQ in some vessels while acetate is never near it, so a rule that is silently analyte-specific would create a spurious group difference. **Do not let a below-LOQ rule differ between the healthy-weight and obesity groups, or between carbohydrate arms.**

4. **Report per-analyte recovery and matrix effects rather than a single method-level figure.** Recent orthogonal comparison shows recovery and sensitivity are compound-specific across SCFAs — no single extraction is optimal for all of them [9]. State recovery per analyte at multiple spike levels and assess matrix effects against the actual fermentation medium, following the matrix-matched validation approach applied across biological matrices in [7]. A single "recovery 85–110%" range collapsed over the panel is not adequate.

5. **Describe the batch structure, randomization, and pooled-QC design explicitly.** State how many batches, how samples were randomized across them, at what step the pooled QC was constituted, how often it was injected, and what RSD thresholds constituted batch acceptance [16]. The constitution step is not a detail: pooling before versus after extraction changes inter-batch correction and can change which analytes appear to differ [17]. Randomization and replication structure are core elements of the community reporting guidance [14]. **For our design specifically, ensure donor and group are not confounded with batch** — all three carbohydrate conditions from a given donor should sit in the same batch wherever possible, so that the donor-paired contrast is within-batch.

6. **Frame the endpoint as a net concentration change in a closed system, and name what that excludes.** State in Methods that the 0→48 h endpoint is net accumulation and that it cannot separate formation from consumption, cross-feeding, or interconversion. Sakata's methodological review is explicit that lumen and fecal concentrations "should not be used as measures of SCFA production or absorption" [10]; the empirical dissociation between fecal concentration and tracer-measured production in humans confirms this [11]; and the *in vivo* absorption and interconversion sinks that a closed vessel lacks are quantified in [12]. Note that this cuts both ways — the closed ex vivo vessel removes host absorption, so our supernatant values are **not** predictive of in vivo luminal or fecal concentrations either.

---

## Prior-claim bullets (mapped to citation numbers)

- **PFBBr derivatization with GC-MS under negative chemical ionization is an established, validated approach for quantifying C1–C5 carboxylic acids in biological matrices, with the derivatization conditions systematically optimized for this analyte class.** [1, 2, 3, 4]
- **Under ECNICI, PFBBr esters of carboxylic acids yield the carboxylate anion as the dominant ion, which is the basis for selected-ion-monitoring quantification and for the method's sensitivity.** [2, 3]
- **Multiple derivatization and derivatization-free GC-MS approaches (PFBBr, propyl chloroformate, BSTFA silylation, direct injection) are in routine use for SCFA quantification in microbiome studies; each has published validation data, and method choice is a documented trade-off rather than a settled question.** [5, 6, 7, 8, 9]
- **Extraction recovery and sensitivity for SCFAs are analyte-specific; no single sample-preparation chemistry is optimal across the whole SCFA panel.** [9]
- **Quantification using stable-isotope-labelled internal standards combined with external calibration curves, with explicit matrix-effect testing, is a validated design for SCFA measurement across biological matrices.** [7, 8]
- **Luminal and fecal SCFA concentrations do not measure the rate of SCFA production and should not be reported as if they do.** [10, 11]
- **In humans, fecal SCFA concentrations dissociate from tracer-measured intestinal SCFA production; plasma concentrations track production better than fecal concentrations do.** [11]
- **A large and analyte-dependent fraction of colonic SCFA is absorbed or interconverted in vivo (systemic availability ≈36% acetate, 9% propionate, 2% butyrate; ~24% acetate→butyrate interconversion), so in vivo and closed-vessel ex vivo measurements are not interchangeable.** [12]
- **Community reporting standards for mass-spectrometry metabolomics require explicit reporting of sample preparation, randomization and replication, calibration, recovery, matrix/ion-suppression assessment, and LOD/LOQ.** [14, 15]
- **Pooled quality-control samples interspersed throughout each analytical batch are the standard mechanism for conditioning, drift correction and precision estimation; how the QC pool is constituted materially affects inter-batch correction and feature selection.** [16, 17]
- **Succinate is the principal precursor of propionate from hexose sugars in the human colonic microbiota, which is why it is measured alongside the three canonical SCFAs in carbohydrate-fermentation panels.** [18, 19]
- **In in vitro fecal batch fermentations, succinate is interpretable as a cross-feeding intermediate between succinate-producing and succinate-utilizing taxa, and its accumulation reflects the balance of that cross-feeding.** [19]
- **Succinate, lactate and formate are routinely reported alongside acetate/propionate/butyrate as fermentation intermediates in gut organic-acid profiling.** [13]
- **Amino-acid-derived (proteolytic) fermentation products are an established analyte class reported alongside saccharolytic products; 5-aminovalerate specifically has been quantified and reported as a gut microbial metabolite, including in panels that also report succinate and the canonical SCFAs.** [20, 21, 22, 23]

---

## Conflicts / nuance

- **The "concentration ≠ production" literature is about in vivo lumen and feces, not closed ex vivo vessels.** Sakata's argument [10] and the tracer studies [11, 12] rest on absorption, viscosity, biofilm and transit — sinks that are largely absent from a sealed batch incubation. Our closed system is therefore *more* defensible than fecal sampling: the 0→48 h change genuinely is net accumulation within a bounded compartment. But it is still not production, because formation, consumption, cross-feeding and interconversion all occur inside the vessel and cannot be separated without isotopic labelling [12, 19]. **The honest framing is: "net accumulation in a closed ex vivo system," not "production," and not "production capacity."** Do not import the fecal critique wholesale as if it invalidated our endpoint, and do not dismiss it as inapplicable either.
- **Sakata explicitly warns against over-emphasizing butyrate** [10]. This is in direct tension with any framing that foregrounds butyrate above acetate and propionate. Given that our claim lock already treats the SDC-versus-RDC butyrate contrast as a primary hypothesis, this warning should be acknowledged in the Discussion as a limitation of emphasis, not buried.
- **PFBBr precedent for SCFAs specifically is thinner than for carboxylic acids generally.** Citation [1] is the strongest direct SCFA precedent, and it is framed as an isotopomer-enrichment assay rather than a concentration assay. Citations [2, 3, 4] establish the chemistry and validation framework for the analyte class but are not SCFA papers. This is a real gap: the Methods should present the workflow as PFBBr/NCI chemistry [2, 3, 4] applied to SCFAs following [1], rather than implying a single canonical SCFA-PFBBr reference method exists.
- **Method-comparison papers disagree about which extraction is best, and say so** [9]. Do not cite any one alternative method paper as establishing a performance benchmark our method must beat; [9] shows the comparison is analyte-conditional. Report our own per-analyte validation figures rather than arguing by reference.
- **5-aminovalerate precedent is precedent for *measurement*, not for a specific biological interpretation.** Citations [21, 22, 23] show it is a legitimately reported analyte in gut-metabolite panels; [20] gives the general proteolytic-fermentation rationale. None establishes a validated interpretation of 5-aminovalerate as a carbohydrate-availability index. Report it descriptively.
- **Citation [22] is the tightest analyte-list precedent but is an infant HMO study from an industry lab.** It is peer-reviewed and appropriate to cite for the panel composition, but should not be leaned on for interpretation.
- **The strongest fecal-vs-production evidence [11] is an inulin study in adults**, not adolescents, and measures in vivo production. It supports the general language guardrail; it does not speak to our population or design.

---

## What we must not overclaim (given claim lock)

Everything in this report is analytical-chemistry and language scaffolding. None of it licenses a biological conclusion. Specifically:

- **Nothing here supports SDC superiority over RDC for any SCFA.** Method citations establish that we can measure the analytes reliably; they say nothing about which carbohydrate yields more of them. Citation [19] shows that substrate chemistry *can* shift the acetate/propionate/butyrate balance in vitro — this is a rationale for asking the question, never evidence for our answer.
- **Nothing here supports obesity-group equivalence or "preserved fermentation capacity."** A validated assay does not convert a null result into demonstrated equivalence. Do not let "our method had good precision and low LOQ" drift into "therefore the absence of a group difference is meaningful." If equivalence is to be claimed at all it requires an equivalence-testing framework, which is outside this prompt's scope.
- **Do not write "production," "produced," "production capacity," "producing capacity," or "SCFA output rate" for our endpoints.** [10, 11] Approved terms: *concentration*, *net accumulation*, *net change (0→48 h)*, *net output*, *accumulated*. Where a rate is unavoidable, write *net accumulation rate* and define it as Δconcentration/Δtime, not as formation.
- **Do not describe succinate as evidence of a specific propionate-pathway flux in our data.** [18] establishes the succinate pathway as the dominant hexose→propionate route across the human microbiota in general. It does not license inferring pathway flux from a single endpoint succinate concentration in our vessels. Succinate accumulation is consistent with, but not diagnostic of, incomplete succinate→propionate conversion.
- **Do not describe 5-aminovalerate as a validated marker of proteolytic fermentation in this study.** [20] supports the general framing that amino-acid-derived metabolites index proteolytic activity; [21–23] support it being a reportable analyte. Report it as an exploratory related fermentation metabolite.
- **Do not use any citation here to support clinical, glycemic, or obesity-treatment relevance.** Citations [11] and [12] are in vivo human physiology studies; they are cited here solely to constrain our measurement language, and must not be repurposed to imply that our ex vivo endpoints carry in vivo physiological meaning.
- **Do not claim our method is superior to alternatives.** [9] specifically undercuts that claim structure. Justify the choice; do not rank it.

---

## Downstream synthesis blocks

### Block 1 — Methods: derivatization and detection

**SECTION:** Methods §2.4 (analytical chemistry)
**PRIOR CLAIM:** PFBBr derivatization with GC-MS under negative chemical ionization is an established, validated approach for quantifying short-chain carboxylic acids in biological matrices.
**CITATIONS (max 3 lead + 2 backup):** Lead — [1] Tomcik 2011 (PMID 21112315); [2] Jia 2003 (PMID 14632119); [3] Tsikas 2016 (PMID 27343144). Backup — [4] Gracia-Moreno 2015 (PMID 25601317); [6] Zhang 2019 (PMID 30683360).
**CONFLICTS / NUANCE:** The direct SCFA-PFBBr precedent [1] is framed as an isotopomer-enrichment assay; the concentration-assay validation framework comes from the broader carboxylic-acid literature [2, 3]. Present the workflow as that chemistry applied to SCFAs, not as a single pre-existing canonical method.
**MANUSCRIPT-SAFE SENTENCE:** "Acetate, propionate, butyrate, succinate and 5-aminovalerate were quantified by gas chromatography–mass spectrometry following pentafluorobenzyl bromide derivatization, with detection under electron-capture negative-ion chemical ionization by selected-ion monitoring of the corresponding carboxylate anion, following established procedures for short-chain carboxylic acids."
**CLAIM-LOCK CHECK:** Supported. Descriptive methods statement; asserts no biological result.

### Block 2 — Methods: calibration, internal standards and QC

**SECTION:** Methods §2.4 (quantification and quality control)
**PRIOR CLAIM:** Stable-isotope internal standards combined with external calibration and explicit matrix-effect testing, with pooled QC samples interspersed across randomized analytical batches, constitute the accepted quantification and QC design for targeted SCFA measurement.
**CITATIONS:** Lead — [7] Rohde 2022 (PMID 35208244); [14] Alseekh 2021 (PMID 34239102); [16] Lewis 2016 (PMID 27479709). Backup — [15] Viant 2019 (PMID 31292445); [17] Ramos 2025 (PMID 39557686).
**CONFLICTS / NUANCE:** How the pooled QC is constituted (pre- vs post-extraction) changes inter-batch behavior [17], so the constitution step must be stated, not assumed. Ensure donor and study group are not confounded with analytical batch.
**MANUSCRIPT-SAFE SENTENCE:** "Quantification used stable-isotope-labelled internal standards added before extraction together with matrix-matched external calibration curves; samples were randomized across analytical batches, and pooled quality-control samples constituted before extraction were injected at regular intervals within each batch to monitor drift and quantify measurement precision, in accordance with community best-practice guidelines for mass-spectrometry-based metabolomics."
**CLAIM-LOCK CHECK:** Supported. Methods description only.

### Block 3 — Methods: units, dilution and below-LOQ handling

**SECTION:** Methods §2.4 (data reporting)
**PRIOR CLAIM:** Reporting standards require explicit statement of units, dilution correction, per-analyte LOD/LOQ, and the rule used for values below the limit of quantification.
**CITATIONS:** Lead — [14] Alseekh 2021 (PMID 34239102); [15] Viant 2019 (PMID 31292445). Backup — [9] Jiang 2025 (PMID 40527045); [8] García-Villalba 2012 (PMID 22865755).
**CONFLICTS / NUANCE:** The below-LOQ rule must be applied identically across analytes and across both study groups and all three carbohydrate arms; an analyte- or group-varying rule can manufacture an apparent difference. This is a live risk for 5-aminovalerate and succinate, which may sit near the LOQ while acetate never does.
**MANUSCRIPT-SAFE SENTENCE:** "Concentrations are reported in µmol/L of fermentation supernatant and corrected for sample dilution; limits of detection and quantification are reported per analyte, and values below the limit of quantification were handled by a single pre-specified rule applied uniformly across all analytes, donors, and carbohydrate conditions, with the number of affected observations reported per analyte and group."
**CLAIM-LOCK CHECK:** Supported. Transparency statement; asserts no result.

### Block 4 — Results/Discussion: endpoint language

**SECTION:** Results and Discussion (SCFA endpoint framing)
**PRIOR CLAIM:** SCFA concentrations measured at a single timepoint, or as a change across an incubation, reflect the net balance of formation and consumption and do not measure production rate.
**CITATIONS:** Lead — [10] Sakata 2019 (PMID 30370625); [11] Kirschner 2025 (PMID 40274191); [12] Boets 2017 (PMID 27510655).
**CONFLICTS / NUANCE:** The cited critique targets in vivo fecal and luminal sampling, where absorption and transit dominate; our closed ex vivo vessel removes those sinks, making net accumulation a cleaner quantity than a fecal concentration. It remains a net quantity, not a production rate, because cross-feeding and consumption continue inside the vessel. State this symmetry rather than citing [10–12] as if they straightforwardly applied.
**MANUSCRIPT-SAFE SENTENCE:** "Because the incubation is a closed system, endpoints are reported as the net change in supernatant concentration between 0 and 48 h; this quantity reflects the balance of microbial formation, consumption and interconversion within the vessel and is not a measure of production rate."
**CLAIM-LOCK CHECK:** Supported, and *protective* of the claim lock — this is the sentence that prevents concentration results from being read as production results.

### Block 5 — Methods/Results: rationale for succinate

**SECTION:** Methods §2.4 (analyte selection) and Results (succinate reporting)
**PRIOR CLAIM:** Succinate is the principal precursor of propionate from hexose sugars in the human colonic microbiota and is routinely reported as a fermentation intermediate in carbohydrate-fermentation panels.
**CITATIONS:** Lead — [18] Reichardt 2014 (PMID 24553467); [19] Zhang 2024 (PMID 39393609). Backup — [13] Tsukuda 2021 (PMID 33723382); [22] Dogra 2021 (PMID 34576834).
**CONFLICTS / NUANCE:** Pathway dominance across the microbiota in general does not license inferring pathway flux from an endpoint succinate concentration in our vessels. Succinate accumulation is consistent with, but not diagnostic of, incomplete succinate→propionate conversion.
**MANUSCRIPT-SAFE SENTENCE:** "Succinate was quantified alongside the three canonical short-chain fatty acids because it is the principal precursor of propionate from hexose sugars in the human colonic microbiota, and its net accumulation is therefore informative about the fate of fermented carbohydrate."
**CLAIM-LOCK CHECK:** Supported as a rationale for measurement. Any statement inferring propionate-pathway flux from our succinate values would be **unsupported**.

### Block 6 — Methods/Results: rationale for 5-aminovalerate

**SECTION:** Methods §2.4 (analyte selection) and Results (5-aminovalerate reporting)
**PRIOR CLAIM:** Amino-acid-derived fermentation products are an established analyte class reported alongside saccharolytic products, and 5-aminovalerate specifically has been quantified as a gut microbial metabolite in panels also reporting succinate and the canonical SCFAs.
**CITATIONS:** Lead — [20] Nyangale 2012 (PMID 23116228); [22] Dogra 2021 (PMID 34576834); [21] Caspani 2022 (PMID 35328680). Backup — [23] Pontifex 2024 (PMID 38835220).
**CONFLICTS / NUANCE:** These establish 5-aminovalerate as a legitimately reported analyte and give the general proteolytic-fermentation rationale. None validates it as an index of carbohydrate availability or of saccharolytic-to-proteolytic switching in a specific dataset.
**MANUSCRIPT-SAFE SENTENCE:** "5-Aminovalerate, an amino-acid-derived microbial fermentation product, was quantified as an exploratory related metabolite alongside the short-chain fatty acids and succinate."
**CLAIM-LOCK CHECK:** Exploratory. Framing 5-aminovalerate as a validated marker of proteolytic fermentation, or interpreting between-group differences in it mechanistically, would be **unsupported**.

### Block 7 — Discussion: limitation on butyrate emphasis

**SECTION:** Discussion (limitations)
**PRIOR CLAIM:** Methodological guidance in the SCFA literature warns against disproportionate emphasis on butyrate relative to the other short-chain fatty acids.
**CITATIONS:** Lead — [10] Sakata 2019 (PMID 30370625). Backup — [13] Tsukuda 2021 (PMID 33723382).
**CONFLICTS / NUANCE:** This sits in direct tension with treating the SDC-versus-RDC butyrate contrast as a primary hypothesis. It should be acknowledged explicitly rather than omitted, and it argues for reporting acetate and propionate with equal prominence in figures and tables.
**MANUSCRIPT-SAFE SENTENCE:** "We report acetate, propionate and butyrate with equal prominence, in line with methodological guidance cautioning against disproportionate emphasis on butyrate in short-chain fatty acid research."
**CLAIM-LOCK CHECK:** Supported. A limitation and reporting-balance statement; strengthens rather than extends the claim set.
