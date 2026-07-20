# Prompt 09 — Statistical frameworks for microbiome differential abundance (ANCOM-BC/BC2, MaAsLin2/3, ALDEx2/3)

**Manuscript target:** Methods 2.8 · Results 3.7 · Discussion 4.6 · **Citation job:** Justify an exploratory, pre-specified single-primary-tool DA strategy for a donor-nested, single-timepoint carbohydrate contrast; explain estimand differences, scale/compositionality handling, and why tools disagree; supply reporting language for relative-vs-qPCR-scaled discordance · **Date:** 2026-07-19

**Verification note.** Every DOI/PMID below was retrieved from PubMed (E-utilities via the bio-research PubMed MCP server) or from Crossref. Nothing here is reconstructed from memory. Publication-status flags for ALDEx3 and MaAsLin 3 are stated explicitly in the table and in §"Version and publication status".

---

## Methods-rationale paragraph (drop-in draft for Methods 2.8)

Differential abundance (DA) testing of 16S amplicon data is not a single well-posed problem: the available frameworks estimate different quantities under different identifying assumptions, and therefore are expected to disagree even when all are correctly implemented [1,2,4,10,11]. Sequence counts are compositional — the library size is set by the instrument, not by the biological system — so no method can recover absolute per-taxon abundance from read counts alone without an explicit or implicit assumption about the unmeasured total microbial load ("scale") [10,15,17]. ANCOM-BC/ANCOM-BC2 target an absolute-abundance log-linear model and address this by estimating a per-sample "sampling fraction" offset and correcting the resulting bias, with ANCOM-BC2 adding multi-group contrasts, covariate adjustment, repeated-measures random effects, and a sensitivity analysis for pseudo-count choice in the presence of structural zeros [1,2]. MaAsLin 2 and MaAsLin 3 target regression coefficients on (typically total-sum-scaled or log-transformed) relative abundance within a generalized linear/mixed model framework; MaAsLin 3 additionally separates the abundance estimand from a distinct prevalence (presence/absence) estimand and can newly ingest an external absolute-scale measurement such as total-load qPCR or spike-ins [3,4]. ALDEx2 targets the centered log-ratio (CLR) difference and propagates technical (Dirichlet–multinomial) sampling uncertainty by Monte Carlo, so its estimand is explicitly ratio-based rather than absolute [5]; recent work reframes every normalization — including the CLR — as an implicit and usually overconfident scale model, and shows that replacing it with an explicit scale model with quantified uncertainty markedly reduces both false positives and false negatives [10,11,16]. The same scale-uncertainty machinery has been extended to nested/repeated-measures designs as scale-reliant mixed-effects models (SR-MEM), implemented in the ALDEx3 R package [12]. Because the present analysis is an exploratory, hypothesis-generating secondary endpoint on a donor-nested, single-timepoint head-to-head carbohydrate contrast (each donor contributing paired RDC and SDC incubations), we pre-specified a single primary DA model with a donor random effect, reported effect sizes with confidence intervals alongside Benjamini–Hochberg-adjusted q values at an exploratory threshold, and ran the remaining frameworks only as pre-specified robustness checks whose purpose is to characterize the sensitivity of each finding to its scale assumption — not to accumulate evidence by counting how many tools agree [7,8,9,19,20,21].

---

## Prioritized citation table

Publication-status column is explicit for every methods paper, per the prompt's requirement.

| # | Title | Year | Journal | DOI/PMID | Status | Why cite |
|---|---|---|---|---|---|---|
| 1 | Analysis of compositions of microbiomes with bias correction (**ANCOM-BC**) | 2020 | Nature Communications 11:3514 | DOI 10.1038/s41467-020-17041-7 · PMID 32665548 | **Peer-reviewed (published)** | Primary ANCOM-BC methods paper. Defines "sampling fraction," models absolute abundance in a linear-regression framework, gives CIs and FDR control. The estimand is absolute-abundance log-fold change, identified via an estimated per-sample offset. |
| 2 | Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures (**ANCOM-BC2**) | 2024 (online 2023-12-29) | Nature Methods 21(1):83–91 | DOI 10.1038/s41592-023-02092-7 · PMID 38158428 | **Peer-reviewed (published)** | Primary ANCOM-BC2 paper — the version to cite for a repeated-measures / donor-nested design. Adds multi-group and ordered-group contrasts, covariate adjustment, random effects, and the pseudo-count sensitivity analysis. (A Research Square preprint of the same work exists, PMID 37205444, DOI 10.21203/rs.3.rs-2778207/v1 — cite the Nature Methods version.) |
| 3 | Multivariable association discovery in population-scale meta-omics studies (**MaAsLin 2**) | 2021 | PLoS Computational Biology 17(11):e1009442 | DOI 10.1371/journal.pcbi.1009442 · PMID 34784344 | **Peer-reviewed (published)** | Primary MaAsLin 2 paper. GLM/GLMM on relative abundances; explicitly evaluated under repeated measures and multiple covariates. Estimand = coefficient on transformed relative abundance. |
| 4 | MaAsLin 3: refining and extending generalized multivariable linear models for meta-omic association discovery | 2026 | Nature Methods 23(3):554–564 | DOI 10.1038/s41592-025-02923-9 · PMID 41540124 | **Peer-reviewed (published) — verified** | Primary MaAsLin 3 paper. **This is now formally published**; cite the Nature Methods version, not the preprint. Critical for this manuscript: MaAsLin 3 separates *abundance* from *prevalence* estimands and accepts external absolute-scale input (qPCR total load or spike-ins) — the exact relative-vs-absolute comparison our Results 3.7 makes. |
| 4b | MaAsLin 3 … (bioRxiv preprint of #4) | 2024 | bioRxiv | DOI 10.1101/2024.12.13.628459 · PMID 39713460 | **Preprint (superseded by #4)** | Listed only so the drafting team does not cite the preprint by mistake. Do not cite. |
| 5 | Unifying the analysis of high-throughput sequencing datasets … by compositional data analysis (**ALDEx2**) | 2014 | Microbiome 2:15 | DOI 10.1186/2049-2618-2-15 · PMID 24910773 | **Peer-reviewed (published)** | Primary ALDEx2 methods paper. Dirichlet Monte Carlo instances + CLR; estimand is a log-ratio difference, not an absolute difference. |
| 6 | Analysis of composition of microbiomes: a novel method for studying microbial composition (**ANCOM**) | 2015 | Microbial Ecology in Health and Disease 26:27663 | DOI 10.3402/mehd.v26.27663 · PMID 26028277 | **Peer-reviewed (published)** | Original ANCOM; needed to explain what ANCOM-BC "corrects." Log-ratio, distribution-free, W-statistic — no p-values or CIs, which is precisely the gap ANCOM-BC filled. |
| 7 | Microbiome differential abundance methods produce different results across 38 datasets | 2022 | Nature Communications 13:342 | DOI 10.1038/s41467-022-28034-z · PMID 35039521 | **Peer-reviewed (published)** | The canonical incomplete-concordance benchmark: 14 methods × 38 datasets, "drastically different numbers and sets of significant ASVs," results depend on preprocessing. ALDEx2 and ANCOM-II most conservative/consistent. (Author Correction: DOI 10.1038/s41467-022-28401-w, PMID 35115546 — cite alongside.) |
| 8 | A realistic benchmark for differential abundance testing and confounder adjustment in human microbiome studies | 2024 | Genome Biology 25:247 | DOI 10.1186/s13059-024-03390-9 · PMID 39322959 | **Peer-reviewed (published)** | Best current benchmark: 19 methods, signals implanted into real profiles, and explicit confounder simulation. Finds only a minority of methods properly control false discoveries; unadjusted confounding produces spurious associations. Strongest support for pre-specifying one adjusted model. |
| 9 | A broken promise: microbiome differential abundance methods do not control the false discovery rate | 2019 | Briefings in Bioinformatics 20(1):210–221 | DOI 10.1093/bib/bbx104 · PMID 28968702 | **Peer-reviewed (published)** | Shows an "alarming excess of false discoveries" across popular DA methods, and that inter-taxon correlation degrades performance. Justifies treating nominal q values as approximate and framing findings as exploratory. |
| 10 | Benchmarking differential abundance analysis methods for correlated microbiome sequencing data | 2023 | Briefings in Bioinformatics 24(1):bbac607 | DOI 10.1093/bib/bbac607 · PMID 36617187 | **Peer-reviewed (published)** | The *correlated-data* benchmark — directly on point for our donor-nested repeated-measures design. "Different DAA-c tools could sometimes produce quite discordant results." Linear-model methods (LinDA, MaAsLin2, LDM) more robust than GLM-based ones; LinDA most robust under strong compositional effects. |
| 11 | LinDA: linear models for differential abundance analysis of microbiome compositional data | 2022 | Genome Biology 23:95 | DOI 10.1186/s13059-022-02655-5 · PMID 35421994 | **Peer-reviewed (published)** | Reference method for the linear-model-on-CLR-with-bias-correction estimand; extends to mixed-effects models for correlated data. Useful as a fourth robustness check and to explain why estimand choice, not implementation quality, drives disagreement. |
| 12 | Scale reliant mixed effects models enhance microbiome data analysis (**SR-MEM; introduces the ALDEx3 R package**) | 2026 | Microbiome 14(1) | DOI 10.1186/s40168-026-02377-x · PMID 41882807 | **Peer-reviewed (published) — verified; ALDEx3 itself has no separate standalone methods paper** | **The correct and only peer-reviewed citation for ALDEx3.** ALDEx3 is a CRAN R package; its methodological content is published here as SR-MEM. Directly on point: mixed-effects modelling with scale treated as a latent variable, able to ingest external qPCR/flow scale measurements, evaluated for longitudinal/hierarchical designs. |
| 13 | Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization | 2025 | Genome Biology 26:139 | DOI 10.1186/s13059-025-03609-3 · PMID 40405262 | **Peer-reviewed (published)** | The primary "scale models" paper. Shows every normalization is an implicit scale assumption; explicit scale models drastically reduce FP and FN. Core citation for why relative-scale and qPCR-scaled models can legitimately disagree. |
| 13b | Beyond Normalization: Incorporating Scale Uncertainty … (preprint of #13) | 2024 | bioRxiv | DOI 10.1101/2024.04.01.587602 · PMID 38617212 | **Preprint (superseded by #13)** | Do not cite; listed to prevent mis-citation. |
| 14 | Explicit Scale Simulation for analysis of RNA-sequencing count data with ALDEx2 | 2025 | NAR Genomics and Bioinformatics 7(3):lqaf108 | DOI 10.1093/nargab/lqaf108 · PMID 40837840 | **Peer-reviewed (published)** | Practical guidance on building scale models in ALDEx2, and the explicit argument that scale models "replace the need for dual cutoff approaches" (arbitrary effect-size + p-value gates). |
| 15 | Microbiome Datasets Are Compositional: And This Is Not Optional | 2017 | Frontiers in Microbiology 8:2224 | DOI 10.3389/fmicb.2017.02224 · PMID 29187837 | **Peer-reviewed (published)** | Standard citation for the compositionality constraint and for why an apparent increase in one taxon's relative abundance can arise from decreases elsewhere. |
| 16 | Incorporating Scale Uncertainty into Differential Expression Analyses Using ALDEx2 | 2026 | Current Protocols 6(2):e70307 | DOI 10.1002/cpz1.70307 · PMID 41637142 | **Peer-reviewed (published)** | Protocol-level, citable how-to for the scale-model workflow; useful in Methods for reproducibility. |
| 17 | Quantitative microbiome profiling links gut community variation to microbial load | 2017 | Nature 551:507–511 | DOI 10.1038/nature24460 · PMID 29143816 | **Peer-reviewed (published)** | Empirical demonstration that relative and load-scaled (quantitative) profiling can yield opposite conclusions — the Bacteroides/Prevotella trade-off as an artefact of relative analysis. Foundational for our relative-vs-qPCR reporting language. |
| 18 | A quantitative sequencing framework for absolute abundance measurements of mucosal and lumenal microbial communities | 2020 | Nature Communications 11:2590 | DOI 10.1038/s41467-020-16224-6 · PMID 32444602 | **Peer-reviewed (published)** | dPCR + 16S absolute-abundance framework; shows absolute (but not relative) measurement revealed a diet-driven decrease in total load and changed per-taxon conclusions. Directly parallels our qPCR-scaled secondary model. |
| 19 | Controlling the false discovery rate: a practical and powerful approach to multiple testing | 1995 | J. R. Stat. Soc. B 57(1):289–300 | DOI 10.1111/j.2517-6161.1995.tb02031.x | **Peer-reviewed (published)** (Crossref-verified; not PubMed-indexed) | The BH procedure itself — cite for the adjustment actually applied. |
| 20 | Statistical significance for genomewide studies (q-values) | 2003 | PNAS 100(16):9440–9445 | DOI 10.1073/pnas.1530509100 · PMID 12883005 | **Peer-reviewed (published)** | Defines the q-value and frames FDR thresholds as a tunable balance of true/false positives rather than a fixed bright line — the standard justification for a more permissive exploratory threshold. |
| 21 | Adjusting for multiple testing — when and how? | 2001 | J. Clin. Epidemiol. 54(4):343–349 | DOI 10.1016/S0895-4356(00)00314-0 · PMID 11297884 | **Peer-reviewed (published)** | The cleanest citable statement that multiplicity adjustment requirements differ between **confirmatory** and **exploratory** analyses, and that the error rate being controlled must be stated. Anchor citation for "exploratory q<0.1 vs confirmatory q<0.05." |
| 22 | The ASA Statement on p-Values: Context, Process, and Purpose | 2016 | The American Statistician 70(2):129–133 | DOI 10.1080/00031305.2016.1154108 | **Peer-reviewed (published)** (Crossref-verified) | Explicitly warns against treating a threshold crossing as a scientific conclusion, and against selective reporting based on which analysis "worked." Backbone of the anti-vote-counting rules. |
| 23 | Vote-counting methods in research synthesis | 1980 | Psychological Bulletin 88(2):359–369 | DOI 10.1037/0033-2909.88.2.359 | **Peer-reviewed (published)** (Crossref-verified) | The canonical statistical demonstration that counting significant vs non-significant results is a biased, low-power synthesis procedure. The methodological principle transfers directly to "N of 4 tools agreed." |
| 24 | Consistent and correctable bias in metagenomic sequencing experiments | 2019 | eLife 8:e46923 | DOI 10.7554/eLife.46923 · PMID 31502536 | **Peer-reviewed (published)** | Measurement bias is multiplicative and taxon-specific; it partly cancels in ratio-based (CLR) estimands but not in absolute-abundance estimands. Explains a *mechanistic* reason relative-scale and qPCR-scaled models can diverge for one genus. |
| 25 | Analysis of microbial compositions: a review of normalization and differential abundance analysis | 2020 | npj Biofilms and Microbiomes 6:60 | DOI 10.1038/s41522-020-00160-w · PMID 33268781 | **Peer-reviewed (published)** | Compact review mapping DA methods to their assumptions and weaknesses; good single citation when Methods needs one review rather than six primaries. |

**Not cited / flagged.** No item above is unverified. There is **no standalone ALDEx3 methods paper** in PubMed or Crossref as of 2026-07-19; if the manuscript needs to cite the software artifact itself, cite the CRAN package entry *in addition to* [12], and label it as software (`⚠ software-only artifact — no independent peer-reviewed methods paper`).

---

## Version and publication status (the item the prompt asked to verify explicitly)

| Framework | Version to cite | Publication status verified 2026-07-19 |
|---|---|---|
| ANCOM | Mandal 2015 [6] | Published (Microb Ecol Health Dis) |
| ANCOM-BC | Lin & Peddada 2020 [1] | Published (Nat Commun) |
| ANCOM-BC2 | Lin & Peddada 2024 [2] | **Published** (Nat Methods 21:83–91). A Research Square preprint also exists; do not cite it. |
| MaAsLin 2 | Mallick 2021 [3] | Published (PLoS Comput Biol) |
| **MaAsLin 3** | Nickols 2026 [4] | **PUBLISHED — Nature Methods 23(3):554–564, DOI 10.1038/s41592-025-02923-9, PMID 41540124.** It is no longer preprint-only; the widely circulated bioRxiv version (10.1101/2024.12.13.628459) is superseded. |
| ALDEx2 | Fernandes 2014 [5]; scale-model extension Nixon 2025 [13] and Gloor 2025 [14] | Published |
| **ALDEx3** | McGovern & Silverman 2026 [12] | **PUBLISHED as a method, but not under its own name.** ALDEx3 is a CRAN R package with **no standalone ALDEx3 methods paper**. Its statistical content — scale-reliant mixed-effects models (SR-MEM) — is published in *Microbiome* 14(1), DOI 10.1186/s40168-026-02377-x, PMID 41882807, where ALDEx3 is described as "an accessible implementation." Cite that paper plus the software artifact; do **not** cite "ALDEx3 (Silverman et al.)" as if a dedicated methods paper existed. |

---

## Prior-claim bullets (each mapped to citation numbers)

- Sequencing counts are compositional; the total is instrument-imposed, so no DA method recovers absolute abundance from counts alone without a scale assumption. [15, 13, 5]
- ANCOM-BC/BC2 estimate a per-sample sampling-fraction offset and target an **absolute-abundance** log-linear estimand with valid p-values, CIs, and FDR control; BC2 extends this to multi-group contrasts, covariates, and repeated measures with a structural-zero/pseudo-count sensitivity analysis. [1, 2]
- MaAsLin 2 targets **relative-abundance** regression coefficients in a GLM/GLMM and was designed and simulation-tested for repeated measures plus covariates. [3]
- MaAsLin 3 splits the estimand into **abundance** and **prevalence** components and can incorporate an **externally measured scale** (qPCR total load or spike-ins); in IBDMDB, 77% of recovered associations were prevalence rather than abundance associations. [4]
- ALDEx2's estimand is a **CLR (log-ratio) difference** with Dirichlet Monte Carlo propagation of technical sampling error; it does not claim to estimate absolute change. [5]
- All normalizations (including CLR and TSS) are implicit scale models; errors in the implicit assumption inflate both false positives and false negatives, and explicit scale models with quantified uncertainty reduce both. [13, 14, 16]
- Scale-uncertainty modelling has been extended to **nested/repeated-measures designs** (SR-MEM, ALDEx3), which can use qPCR or flow-cytometry load as an external scale input and controlled FDR in simulations and reanalyses. [12]
- DA tools produce **incompletely concordant** results on the same data: drastically different numbers and sets of significant features across 14 methods × 38 datasets, with results also depending on preprocessing. [7]
- On realistic simulations with implanted signals and confounders, only a subset of 19 DA methods properly controlled false discoveries; unadjusted confounding produced spurious associations in real data. [8]
- Many popular DA methods show excess false discoveries; inter-taxon correlation makes this worse. [9]
- For **correlated / repeated-measures** microbiome data specifically, DA tools "could sometimes produce quite discordant results," with linear-model methods more robust than GLM-based ones. [10, 11]
- Relative and load-scaled profiling can yield **opposite** biological conclusions; a well-known taxonomic trade-off was shown to be an artefact of relative analysis. [17, 18]
- Measurement bias is consistent, multiplicative, and taxon-specific — so it behaves differently in ratio-based versus absolute-scale estimands. [24]
- Multiplicity adjustment expectations legitimately differ between **confirmatory** and **exploratory** analyses, and the error rate controlled must be stated explicitly. [21, 19, 20]
- Threshold-crossing is not a scientific conclusion, and reporting should not be conditioned on which analysis produced significance. [22]
- Counting significant versus non-significant results across analyses is a biased and low-power synthesis procedure. [23]

---

## Conflicts / nuance

1. **Nearing et al. [7] recommend a consensus-of-multiple-tools approach; the scale-model literature [13, 12] implies this is not principled.** Nearing's recommendation is a *pragmatic robustness heuristic* offered in the absence of a way to adjudicate between estimands. Nixon et al. [13] show that disagreement between tools is largely disagreement about the unmeasured scale, so intersecting tool outputs does not average away the error — it silently selects for features robust to *all* the implicit scale assumptions simultaneously, with an unknown resulting error rate. **These are not the same recommendation and should not be blended.** Our manuscript should present the intersection descriptively while explicitly declining to treat it as a test.
2. **Benchmarks disagree with each other about which tool wins.** Nearing [7] favours ALDEx2/ANCOM-II for consistency; Yang & Chen [10] favour LinDA/MaAsLin2/LDM for correlated data; Wirbel [8] favours classical linear models, limma, and fastANCOM. There is no benchmark-supported claim that any one tool is universally best — which is itself the argument for pre-specification rather than post-hoc tool selection.
3. **"Absolute" is a modelled quantity, not an observation.** ANCOM-BC's absolute-abundance estimand is identified through an *estimated* sampling fraction, not a measured load [1]. Our qPCR-scaled model uses a *measured* total-load proxy, which is a different (and, in principle, better-identified) route to the same estimand — but qPCR total 16S load is itself subject to 16S copy-number variation and extraction bias [24]. Neither is ground truth.
4. **MaAsLin 3's prevalence estimand can generate apparent "new" findings.** A genus can be significant for prevalence and null for abundance, or vice versa [4]. If we run MaAsLin 3, we must report which estimand a hit belongs to; otherwise the hit count is not comparable to ANCOM-BC2 or ALDEx2 output.
5. **Excess-FDR findings [9] and the realistic-benchmark findings [8] both argue that nominal q values are optimistic.** This cuts against, not for, a permissive threshold. The defensible position is a permissive threshold *paired with* explicit exploratory framing and no confirmatory language — not a permissive threshold used to claim discovery.
6. **n = 16 donors, single timepoint.** Nearing [7] found that the number of features detected correlates with sample size, depth, and effect size. At n = 16 with a within-donor paired contrast, power is limited and the honest expectation is few or zero surviving features; a null DA result is a legitimate and reportable outcome, not a failure of the analysis.
7. **A CLR-based result and a qPCR-scaled result answering different questions is the expected behaviour of correctly functioning software** [13, 15, 17]. Discordance is informative about scale, not evidence that one model is broken.

---

## What we must not overclaim (given claim lock)

- **Do not** use any DA result to support SDC superiority over RDC. A genus-level DA hit on a paired SDC-vs-RDC contrast is a compositional/abundance observation about the community, not evidence of a superior fermentation outcome, and it cannot be used to reinforce an SCFA superiority claim.
- **Do not** claim *absolute* expansion of *Fusicatenibacter* (or any genus) on the basis of a relative-abundance model. Only a qPCR-scaled or otherwise scale-aware model addresses absolute change [1, 4, 12, 17, 18], and if that model is null, the absolute claim is unsupported.
- **Do not** describe any genus as an SDC utilizer, a butyrate producer, or the mechanistic source of a butyrate signal on the basis of DA co-occurrence. DA models estimate association with a covariate; they contain no metabolic information.
- **Do not** frame obesity-group similarity in DA output as "preserved fermentation capacity" or equivalence. Failure to reject is not equivalence, especially at n = 8 per group with acknowledged excess-FDR behaviour [9] and sample-size-dependent detection [7].
- **Do not** state or imply that agreement across two or more DA tools raises confidence in a quantifiable way, or report "significant by k of n methods" as though k were a statistic [23, 22, 13].
- **Do not** report an exploratory q<0.1 hit with confirmatory vocabulary ("demonstrates," "establishes," "confirms," "shows that SDC increases"). Use "nominally associated," "exploratory," "hypothesis-generating."
- **Do not** present a DA finding as a responder-phenotype marker, a predictor of dietary response, or a clinical/glycemic signal. Nothing in this design supports that, and no citation in this table can be used to bridge to it.
- **Do not** silently switch the reported primary model to whichever framework produced significance. Pre-specification must be stated and honoured, and all pre-specified frameworks reported regardless of outcome [22, 8].

---

## Five interpretive rules (for Methods 2.8 / Results 3.7 / Discussion 4.6)

**Rule 1 — Name the estimand before the p-value.** Every DA statement in the manuscript must specify what quantity was estimated: a CLR log-ratio difference (ALDEx2), a relative-abundance regression coefficient (MaAsLin 2), a modelled absolute log-fold change via estimated sampling fraction (ANCOM-BC2), a measured-scale absolute coefficient (MaAsLin 3 with qPCR / ALDEx3 SR-MEM), or a prevalence coefficient (MaAsLin 3). Two frameworks disagreeing about the same genus is not a contradiction if they are estimating different things. [1, 2, 3, 4, 5, 12, 13]

**Rule 2 — One pre-specified primary model; everything else is a labelled sensitivity analysis.** The primary DA model is fixed in advance (donor random effect, carbohydrate condition as the fixed effect of interest, single 48 h timepoint), and its result is *the* result. Additional frameworks are reported in full — including where they disagree — as sensitivity analyses characterizing dependence on the scale assumption. The primary model is never re-chosen after seeing results. [8, 10, 22]

**Rule 3 — Exploratory thresholds are declared, not discovered.** Report BH-adjusted q values, state the threshold in Methods before Results, and use q<0.10 explicitly labelled as an exploratory, hypothesis-generating screen for this secondary endpoint — distinguished from the q<0.05 (or α=0.05) confirmatory standard applied to the primary SCFA endpoints. State which error rate is being controlled. Given documented excess-FDR behaviour in DA methods, treat q values as approximate and never as a discovery certificate. [19, 20, 21, 9, 8]

**Rule 4 — Relative-scale significance without absolute-scale significance is reported as scale-dependence, not as a finding weakened by a failed replication.** Recommended template sentence: *"Genus X showed a higher relative abundance in SDC than RDC incubations (ANCOM-BC2 / MaAsLin 2 log-fold change …, 95% CI …, q = …), but the corresponding qPCR-scaled absolute-abundance model did not detect a difference (log-fold change …, 95% CI …, q = …). Because sequencing counts are compositional, a relative-abundance shift is consistent with either an increase in genus X or a decrease in other community members, and the present data do not distinguish these; we therefore describe this as an exploratory relative-scale association and make no claim of absolute expansion."* Always give the effect size and CI for the *null* absolute-scale model, not just its q value — the width of that interval, not the failure to cross a threshold, is what the reader needs. [1, 4, 12, 13, 15, 17, 18, 24]

**Rule 5 — Never count methods.** Do not write "significant in 3 of 4 tools," do not tabulate a vote tally as evidence, and do not define a "consensus" set and then test it. Tool agreement is descriptive; report the intersection and the disagreements transparently, and interpret disagreement as information about which assumption a finding depends on. If a robustness statement is wanted, phrase it as dependence: *"this association was present in ratio-based models but not in scale-aware models, indicating it is sensitive to the assumption of constant total microbial load."* [23, 22, 13, 7]

---

## Downstream synthesis blocks

**SECTION:** Methods 2.8 (Differential abundance analysis)
**PRIOR CLAIM:** DA frameworks estimate different quantities under different scale assumptions, and for correlated/repeated-measures designs the tools are known to give discordant results; therefore a single pre-specified donor-aware primary model with declared exploratory FDR threshold is the defensible design.
**CITATIONS (max 3 lead + 2 backup):** Lead — [2] ANCOM-BC2 (10.1038/s41592-023-02092-7); [13] Nixon 2025 scale models (10.1186/s13059-025-03609-3); [10] Yang & Chen correlated-data benchmark (10.1093/bib/bbac607). Backup — [4] MaAsLin 3 (10.1038/s41592-025-02923-9); [21] Bender & Lange (10.1016/S0895-4356(00)00314-0).
**CONFLICTS / NUANCE:** Benchmarks do not converge on a single best tool ([7] vs [10] vs [8]); this supports pre-specification rather than tool-shopping. Nominal FDR control is optimistic in practice [9].
**MANUSCRIPT-SAFE SENTENCE:** "Because differential abundance frameworks target different estimands under different assumptions about the unmeasured total microbial load, and because published benchmarks report incomplete concordance among them — particularly for correlated, repeated-measures designs — we pre-specified a single primary model (donor as a random effect, carbohydrate condition as the fixed effect, 48 h only), reported Benjamini–Hochberg q values against a declared exploratory threshold of q<0.10, and analysed the remaining frameworks as labelled sensitivity analyses."
**CLAIM-LOCK CHECK:** Supported (methodological statement only; makes no biological claim).

---

**SECTION:** Methods 2.8 (scale handling / qPCR)
**PRIOR CLAIM:** Relative-abundance and total-load-scaled models answer different questions, and external scale measurement (qPCR) can be incorporated into modern DA frameworks.
**CITATIONS:** Lead — [4] MaAsLin 3 (10.1038/s41592-025-02923-9); [12] McGovern & Silverman SR-MEM/ALDEx3 (10.1186/s40168-026-02377-x); [15] Gloor 2017 (10.3389/fmicb.2017.02224). Backup — [17] Vandeputte 2017 (10.1038/nature24460); [18] Barlow 2020 (10.1038/s41467-020-16224-6).
**CONFLICTS / NUANCE:** qPCR total 16S load is a proxy subject to copy-number and extraction bias [24]; "absolute" in ANCOM-BC is modelled, not measured [1]. Neither route is ground truth.
**MANUSCRIPT-SAFE SENTENCE:** "In parallel with relative-abundance models, we fitted absolute-scale models in which per-sample total bacterial 16S qPCR load was supplied as an external scale estimate, since relative and load-scaled profiling are known to support different conclusions from the same sequencing data."
**CLAIM-LOCK CHECK:** Supported.

---

**SECTION:** Results 3.7 (relative vs qPCR-scaled discordance)
**PRIOR CLAIM:** A genus significant on relative abundance but null on the qPCR-scaled model is a scale-dependent, exploratory observation, not evidence of absolute expansion.
**CITATIONS:** Lead — [15] Gloor 2017; [17] Vandeputte 2017; [13] Nixon 2025. Backup — [4] MaAsLin 3; [18] Barlow 2020.
**CONFLICTS / NUANCE:** The null absolute-scale model at n = 16 may simply be underpowered; report its effect size and CI so the reader can see the interval rather than inferring absence of effect from a q value. Failure to reject is not equivalence.
**MANUSCRIPT-SAFE SENTENCE:** "[Genus] showed a higher relative abundance in SDC than in RDC incubations (log-fold change …, 95% CI …, q = …), whereas the qPCR-scaled absolute-abundance model did not detect a corresponding difference (log-fold change …, 95% CI …). Because amplicon counts are compositional, this pattern is equally consistent with an increase in [genus] and with a decrease in other taxa, and we therefore report it as an exploratory, scale-dependent association rather than evidence of absolute expansion."
**CLAIM-LOCK CHECK:** Exploratory. Becomes **unsupported** if written as absolute expansion, as SDC superiority, or as a mechanistic butyrate attribution.

---

**SECTION:** Results 3.7 (cross-method reporting)
**PRIOR CLAIM:** Cross-framework agreement and disagreement are reported descriptively; they are not combined into a test.
**CITATIONS:** Lead — [7] Nearing 2022 (10.1038/s41467-022-28034-z, + correction 10.1038/s41467-022-28401-w); [23] Hedges & Olkin 1980 (10.1037/0033-2909.88.2.359); [22] Wasserstein & Lazar 2016 (10.1080/00031305.2016.1154108). Backup — [8] Wirbel 2024; [13] Nixon 2025.
**CONFLICTS / NUANCE:** Nearing explicitly recommends a consensus approach; we adopt the descriptive part (report multiple methods transparently) while declining the inferential part (treating agreement as evidence), on the grounds that disagreement reflects differing scale assumptions rather than independent evidence [13].
**MANUSCRIPT-SAFE SENTENCE:** "Results from all pre-specified frameworks are reported in Supplementary Table S…; agreement across methods is presented descriptively and was not used as a criterion for significance, since the frameworks differ in their estimands and scale assumptions rather than constituting independent tests."
**CLAIM-LOCK CHECK:** Supported.

---

**SECTION:** Discussion 4.6 (limitations of the DA analysis)
**PRIOR CLAIM:** With n = 16 donors, a single timepoint, and documented excess-FDR behaviour among DA methods, genus-level results here are hypothesis-generating and require independent confirmation on a scale-aware design.
**CITATIONS:** Lead — [9] Hawinkel 2019 (10.1093/bib/bbx104); [8] Wirbel 2024 (10.1186/s13059-024-03390-9); [7] Nearing 2022. Backup — [10] Yang & Chen 2023; [21] Bender & Lange 2001.
**CONFLICTS / NUANCE:** Detection counts scale with sample size, depth, and effect size [7], so both positive and null genus-level results at this n must be interpreted cautiously; a null is not equivalence and a hit is not a discovery.
**MANUSCRIPT-SAFE SENTENCE:** "Genus-level differential abundance was a secondary, exploratory endpoint. Given the number of donors, the single 48 h timepoint, and evidence that differential abundance methods can exceed their nominal false-discovery rates and detect feature counts that track sample size and sequencing depth, these associations are hypothesis-generating and would require independent confirmation in a design with measured microbial load and larger donor numbers before any mechanistic interpretation."
**CLAIM-LOCK CHECK:** Supported.

---

**SECTION:** Discussion 4.6 (why methods disagree — interpretive framing)
**PRIOR CLAIM:** Disagreement among DA frameworks is expected and interpretable as disagreement about the unmeasured biological scale, not as software error or as a reason to prefer whichever result is significant.
**CITATIONS:** Lead — [13] Nixon 2025; [12] McGovern & Silverman 2026; [15] Gloor 2017. Backup — [24] McLaren 2019; [25] Lin & Peddada 2020 review.
**CONFLICTS / NUANCE:** ALDEx3/SR-MEM is published as SR-MEM in *Microbiome* [12]; ALDEx3 itself has no dedicated methods paper and should be cited as software alongside [12].
**MANUSCRIPT-SAFE SENTENCE:** "The frameworks applied here differ chiefly in what they assume about the unmeasured total microbial load; consequently, where relative-scale and load-scaled models diverge, we interpret the divergence as indicating sensitivity to that assumption rather than as a discrepancy to be resolved in favour of either result."
**CLAIM-LOCK CHECK:** Supported.
