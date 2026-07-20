# Deep-research citation package — synthesis index

**Project:** Abbott–Lurie ZC07, Project 2 · **Generated:** 2026-07-19
**Source prompts:** `docs/deep-research-citation-prompts.2026-07-19.md` (18 prompts, all executed)
**Shared guardrails:** `research/_SHARED_CONTEXT.md`

All 18 reports are complete. Every citation was retrieved from a tool (PubMed E-utilities, Consensus with PubMed confirmation, Crossref, or a publisher record); items that could not be verified are flagged `⚠ UNVERIFIED` in place and are listed in §5 below.

---

## 1. Report map

| Prompt | File | Manuscript target |
|---|---|---|
| 01 | `prompt-01-pediatric-obesity-heterogeneity.md` | Introduction ¶1 |
| 02 | `prompt-02-scfa-net-accumulation.md` | Introduction ¶2; Discussion 4.3 |
| 03 | `prompt-03-carb-structure-digestibility.md` | Introduction ¶3; Discussion 4.2–4.3 |
| 04 | `prompt-04-exvivo-fermentation-assay.md` | Introduction ¶4–5; Discussion 4.7–4.8 |
| 05 | `prompt-05-obesity-fermentation-capacity.md` | Introduction ¶5; Discussion 4.4 |
| 06 | `prompt-06-stool-collection-inoculum.md` | Methods 2.1–2.3 |
| 07 | `prompt-07-scfa-gcms-quantification.md` | Methods 2.4 |
| 08 | `prompt-08-16s-qpcr-scaling.md` | Methods 2.5–2.6; Discussion 4.6, 4.8 |
| 09 | `prompt-09-differential-abundance-methods.md` | Methods 2.8; Results 3.7; Discussion 4.6 |
| 10 | `prompt-10-nested-mixed-models.md` | Methods 2.8; Results 3.5–3.6 |
| 11 | `prompt-11-interpersonal-heterogeneity.md` | Discussion 4.5 (centerpiece) |
| 12 | `prompt-12-background-fermentation-controls.md` | Discussion 4.3 |
| 13 | `prompt-13-composition-function-discordance.md` | Discussion 4.4–4.5 |
| 14 | `prompt-14-fusicatenibacter.md` | Discussion 4.6; Results 3.7–3.8 |
| 15 | `prompt-15-rodent-translation.md` | Discussion 4.2; Introduction ¶3–4 |
| 16 | `prompt-16-mechanistic-followup.md` | Discussion 4.6–4.8; Conclusion |
| 17 | `prompt-17-adolescent-microbiome.md` | Introduction ¶1 & ¶4; limitations |
| 18 | `prompt-18-competitive-prior-art.md` | Introduction ¶3; Discussion 4.2–4.3 |

---

## 2. Findings that change the manuscript's framing

These are results the prompts did not anticipate. Each is load-bearing.

### 2.1 *Fusicatenibacter* is not a described butyrate producer (prompts 14, 16, 02, 03)

The type-strain description (Takada et al. 2013, PMID 23625266) reports fermentation end products from glucose as **lactic, formic, acetic and succinic acid — butyrate is not among them.** No isolate, pure-culture, genomic (`but`/`buk`), or flux study establishes butyrogenic capacity for any strain in this genus. The common "butyrate producer" label propagates from Lachnospiraceae family membership, not data; one traced instance (Wikipedia genus entry) miscites the claim directly to the 2013 paper that contradicts it.

**Consequence:** the manuscript should not merely *hedge* on butyrogenesis — it should positively note the genus is not a described butyrate producer. Vital et al. 2014 (PMID 24757212) adds that butyrate synthesis is polyphyletic and not inferable from 16S. The defensible reframe is a **cross-feeding-donor hypothesis** (lactate/acetate → butyrogenic *Anaerostipes*/*E. hallii*), explicitly marked untested; Belenguer 2006 (PMID 16672507) supplies the ¹³C co-culture design that would test it.

Supporting context, not mechanism: *Fusicatenibacter* is carbohydrate-responsive (resistant maltodextrin roughly doubled it by qPCR, PMID 35683992), and two independent studies show it shifting **without** concurrent SCFA change — so a relative signal lacking qPCR-scaled confirmation is the previously observed norm, not an anomaly. Jin 2019 (PMID 30919578) is the only safe correlative framing. Closest demographic match (adolescents with obesity, PMID 35344283) found it *higher* in NAFLD, reversing the beneficial-taxon valence.

### 2.2 The rodent anchor's butyrate contrast was itself null (prompt 15)

Reading the Plaza-Díaz PDF in full: butyrate was **351.9 / 478.3 / 618.8 µg/g with no significant difference.** Our nonsignificant SDC–RDC butyrate result therefore fails to replicate nothing. Given commit `5d1e94e` elevates butyrate to a primary hypothesis, the stated rationale needs revisiting — prespecification, not the anchor paper, must carry that choice.

The anchor also **confounds digestion rate with fermentable-fiber substitution**: at equal declared total fiber, HFD cellulose is replaced by resistant maltodextrin + inulin:FOS, so its acetate/propionate rise cannot be attributed to digestion rate. It reports *decreased* α-diversity in the treated arm and never mentions *Fusicatenibacter*.

### 2.3 "Slow-digesting" was never a colonic-substrate label (prompts 03, 15)

Two independent derivations of the same conclusion:

- **Definitional (03):** RDS/SDS/RS is defined purely by *host enzymic hydrolysis rate*; Englyst offers only the resistant fraction as a colonic-delivery guide. INFOGEST 2.0 authors state the static protocol is unsuitable for simulating digestion *kinetics*, so no static pre-digestion step can restore the property distinguishing SDC from RDC.
- **Empirical (15):** human ileostomy data (PMID 20211041) show isomaltulose — 26.4% of the ISR blend — is 95–99% digested and 94–96% absorbed in the small intestine. In vivo it is largely a glycemic manipulation that barely reaches the colon; our assay presents it intact.

**Consequence:** a null SDC-vs-RDC contrast in a direct-exposure assay is the *expected* result, not a disappointment. In vivo digestion rate and ex vivo substrate availability are different independent variables.

### 2.4 The local "preprint" PDF is a published, Abbott-affiliated paper (prompt 18)

`research/impact slow v rapdi carb pre pub children randomized cross over.pdf` = **Gillen et al. 2021, Clin Nutr, PMID 34130017, NCT03185884** — published, not a preprint. Same RDC/SDC contrast, pre-pubertal children, but measures only host substrate oxidation in healthy-weight children.

**Novelty this licenses:** with Plaza-Díaz covering rodent microbiome and Gillen covering human host, the **human-microbial cell of that 2×2 is empty.**

### 2.5 Near-exact prior art exists and is in-house (prompts 04, 06, 17, 18)

Holmes et al. 2020 (mBio, PMID 32788375) — ex vivo prebiotic fermentation, 17 donors aged 10–18 with obesity, donor-resolved SCFA + 16S, **P.C. Seed co-author.** Surfaced by four independent agents. PubMed keyword search misses it (says "in vitro system," not "in vitro fermentation"); only Consensus found it.

It has **no healthy-weight arm**, and the only pediatric study with a weight contrast (Qiu 2025) **pooled** donors. Precise licensed claim: *donor-resolved, weight-stratified adolescent fermentation has not been reported.*

Critically, Holmes reports fecal SCFA concentration, ex vivo production capacity, and obesity markers did **not** positively correlate — which makes it unusable for any preserved-capacity framing, and argues for attaching **no health valence** to any observed SCFA difference.

---

## 3. Cross-cutting tensions to resolve before drafting

### 3.1 Functional redundancy vs donor heterogeneity — surfaced independently by 6 agents

The single most important unresolved tension, because donor heterogeneity is the Discussion 4.5 centerpiece.

- **Against heterogeneity:** Reichardt 2017 (PMID 29192904) — SCFA output "surprisingly reproducible" across compositionally divergent donors. Reichardt 2018, Oliver 2021, Yao 2024 concur.
- **For analyte-specific reconciliation:** Fan 2026 — fiber reliably raises total SCFA and acetate, but **propionate and butyrate vary substantially between individuals.** Redundancy dominates for total SCFA and simple substrates; donor-specificity dominates for butyrate and complex/food-matrix substrates.
- **Guard against the easy escape:** Ho & Huang 2025 (Cell Systems) — the statistical signature of functional redundancy can arise from **averaging alone**, so redundancy must never be asserted as demonstrated, and cannot be used to explain our obesity null as "preserved capacity."

**Recommended stance:** frame analyte-specifically. This maps directly onto our own analyte-specific results.

### 3.2 pH is undocumented, not known-absent — **needs a wet-lab answer**

Four reports (03, 04, 07, 12) make pH the load-bearing mechanism for butyrate/propionate stoichiometry (Duncan 2009 PMID 19397676, Walker 2005 PMID 16000778, LaBouyer 2022). In an uncontrolled-pH batch system, "carbohydrate effect on butyrate" and "acid-accumulation effect on butyrate" are partially entangled.

One agent asserted pH was "neither controlled nor measured." **This was checked and is not established.** `manuscript/manuscript-outline.md:217` lists pH among parameters still *to be reported*; `docs/gut-microbes-claim-audit.2026-07-19.md:210` says the repo lacks manuscript-ready documentation for medium, inoculum, atmosphere and pH. If the vessels were buffered, the caveat weakens substantially.

**Action: obtain the actual fermentation pH protocol from the bench team.** Four reports depend on it, and butyrate is now a primary hypothesis.

### 3.3 The concentration-is-not-production caveat cuts both ways (prompts 02, 05, 07)

Sakata 2019 (PMID 30370625) states luminal/fecal SCFA concentration must not be read as production or absorption. Prompt 07 supplies the correction: that critique targets **in vivo** sampling, where absorption and transit are the dominant sinks. A closed ex vivo vessel removes those sinks, making net 0→48 h accumulation a *cleaner* quantity than fecal concentration — still not a rate. State the asymmetry explicitly; do not import the fecal critique wholesale, and do not deploy it only where it flatters the design.

Sakata separately warns against disproportionate butyrate emphasis — in tension with commit `5d1e94e`. Written up as an acknowledged limitation.

### 3.4 Multiple-DA-method use: descriptive yes, inferential no (prompts 08, 09)

Nearing 2022 is the obvious citation for running several DA tools but recommends a **consensus** approach, contradicting the triangulation-without-concordance stance. Resolution adopted: take Nearing's *descriptive* claim (tools disagree), decline the *inferential* one, on the grounds from Nixon 2025 that disagreement reflects unmeasured scale — intersecting outputs does not average away error and has unknown error rate. This is what makes "never count methods" defensible rather than merely convenient.

---

## 4. Claim-lock hazards — citations that look supportive but are forbidden

Each was found and flagged by an agent rather than silently used.

| Citation | Why it is a hazard | Permitted use |
|---|---|---|
| **Li et al. 2015** (PMID 26256258) | Concludes resistant starch "may benefit individuals across BMIs" from similar lean-vs-obese SCFA — the exact forbidden equivalence claim. Conflicts with Aguirre 2014 and Vieira 2021. | Cite only as contrasting evidence; never echo the conclusion. |
| **Bartsch 2025** (J Nutr) | Reports arabinoxylan raising *Fusicatenibacter* at **q = 0.18** — not significant at any threshold, but reads as superficially supportive of the guardrailed claim. | Baseline-community-type finding only. Never for a *Fusicatenibacter* claim. |
| **Gerasimidis 2022** | "Similar fermentative capacity despite dysbiosis" is structurally the forbidden obesity-equivalence claim. | Design precedent only. |
| **Reichardt 2018** (ISME J) | Best ex vivo donor-as-unit design precedent, but its stated conclusion is *functional redundancy*. | Cite for design; never for capacity claims. |
| **Holmes 2020 / 2022** | Emphasize donor-dependent responses and "personalized prebiotics" — one step from locked responder-phenotype and dietary-prediction claims. | Justifies donor-aware variance modeling only. |
| **Gurry 2021** | Prediction result is ex vivo *self*-prediction only. | Flagged non-extensible. |
| **Zeevi 2015** | Precision-nutrition prediction. | Validation *architecture* only; never as endorsement of microbiome-based dietary-response prediction. |

**Version traps** (would produce uncitable references):
- **MaAsLin 3** is published: Nickols et al., *Nat Methods* 23(3):554–564 (2026), PMID 41540124. The bioRxiv version is superseded — do not cite it.
- **ALDEx3** has *no* standalone methods paper. Its statistics are published as scale-reliant mixed-effects models (SR-MEM): McGovern & Silverman, *Microbiome* (2026), PMID 41882807, where ALDEx3 is "an accessible implementation." Cite that paper **plus** the software artifact.
- **ANCOM-BC2** has a Research Square preprint twin (PMID 37205444). Correct reference is *Nat Methods* 21(1):83–91, **PMID 38158428**.

---

## 5. Items flagged `⚠ UNVERIFIED` — resolve before submission

| Report | Item | Status |
|---|---|---|
| 07 | Kirschner 2023, Clin Nutr ESPEN conference abstract | Superseded by peer-reviewed Kirschner 2025 AJCN (PMID 40274191). **Do not cite the abstract.** |
| 10 | Shrout & Fleiss 1979 (ICC) | Not verified this run. Check manually if added. |
| 10 | Stoffel, Nakagawa & Schielzeth 2017, *MEE* (rptR article) | Not verified; distinct from the verified CRAN record. |
| 13 | Kelly et al. 2015, PERMANOVA power, *Bioinformatics* | No identifier retrieved. |
| 14 | *Fusicatenibacter faecihominis* (2nd species) | Name and valid-publication status via LPSN only. |
| 15 | Lina et al. 2002, *Food Chem Toxicol* | Backup only; unused in any synthesis block. |
| 08 | ALDEx2 2014 DOI (PMID 24910773 confirmed) | PMID confirmed, DOI not fetched. |
| 03 | Englyst 1992 (PMID 1330528) | No DOI indexed in PubMed — cite by PMID. This is expected, not an error. |

---

## 6. Statistical and methodological cautions

- **ICC robustness does not transfer.** Schielzeth 2020's LMM-robustness result applies to *fixed-effect estimates*, not to the random-effect variance components and ICCs reported in §3.5. This is the specific justification for bootstrap CIs on variance components.
- **Dispersion check required.** Anderson & Walsh 2013: the modest obesity-group PERMANOVA term (R² ≈ 0.026–0.055) must be paired with a betadisper check before it can be read as a centroid shift. Confirm this was run at 48 h.
- **Below-LOQ rule must be uniform.** Succinate and 5-aminovalerate may sit near LOQ while acetate never does; an analyte-varying rule could manufacture an apparent group difference.
- **Freeze–thaw splits by endpoint.** Mild for DNA-based composition below ~4 cycles, but demonstrably alters quantitative SCFA output from cryopreserved inocula (Aguirre 2015). Belongs in Methods *and* Limitations; "acceptable alternative to fresh" must never become "equivalent."
- **Endpoint may under-resolve rate.** Kaur 2011, Rose 2009, Gu 2018 all locate slow-vs-fast differences in the *shape* of the 0–48 h curve. A single 48 h net-accumulation endpoint is a resolution limit — this both belongs in limitations and pre-empts over-reading the null.
- **qPCR scaling is a partial load correction.** With n=16, the absolute-scale null is absence of evidence, not evidence of absence — and the mirror-image overclaim (proven non-expansion) is equally forbidden.
- **No consensus reporting checklist exists** for anaerobic batch fecal fermentation. Recommended composite: Pérez-Burillo 2021 (fermentation parameters) + STORMS + MIxS/MIMARKS.
- **Terminology:** "estimated 16S gene-copy equivalents," not cell counts (Louca 2018 on unresolved copy-number correction; Johnson 2019 on intragenomic variation; Props 2017 — relative enrichment does not imply absolute outgrowth).

---

## 7. Strongest single citations by job

| Job | Citation |
|---|---|
| Net accumulation ≠ production | Sakata 2019, PMID 30370625 |
| SCFA vessel ≠ host exposure | Boets 2017 — 36% / 9% / 2% systemic availability |
| Composition → function fails empirically | Boets 2017 — colonic acetate→butyrate conversion unrelated to fecal butyrogenic capacity |
| Batch-fermentation SOP | Pérez-Burillo 2021, PMID 34089022 |
| Background/control-well accumulation | Van den Abbeele 2022, PMID 35456812 — ~20 mM acetate, ~10 mM butyrate with no added carbohydrate |
| Fresh vs cryopreserved inoculum | Aguirre 2015, J Microbiol Methods |
| Donor as experimental unit | Lazic 2018 (pseudoreplication); Lazic 2010 (replicate averaging) |
| Median-split instability | Hecksteden 2018, PMID 29357481 — 11/20 consistently classified; near-exact analogue of our 7–12/16 |
| Prediction from concurrent features fails | Sze & Topçuoğlu 2019, PMID 31266879 — RF explained ≤14% of fecal SCFA variance |
| Overfitting warrant | Varma & Simon 2006, PMID 16504092 |
| Gold-standard mechanism stack we lack | Ze et al. 2012 (*R. bromii* keystone), PMID 22343308 |
| Adolescent heterogeneity anchor | Ryder 2019 — BMI change −50.2% to +12.9% |
| Analyte-specific SCFA under varied starch | Giuberti 2013; Deehan 2020 |
| Background medium masks substrate effect | Poppe 2023 |
| Propionate pathways often unstimulated | Van-Wehle & Vital 2024, PMID 39080292 |

---

## 8. Genuine gaps found (defensible novelty, and its limits)

1. **Pediatric ex vivo fermentation is near-empty.** `adolescent[TIAB]` + in vitro/ex vivo fermentation returns **0** PubMed hits; MeSH intersection returns 6; adolescent-obesity *compositional* studies return 26. Caveat: keyword search missed Holmes 2020, so raw counts modestly understate the corpus.
2. **No published comparator for cross-substrate rank correlation.** No prior study reports a formal cross-substrate rank correlation of SCFA response with a CI, so ρ ≈ 0.83–0.96 is novel but externally unanchored.
3. **Human-microbial RDC/SDC cell is empty** (see §2.4).
4. **What the gaps do *not* license:** the assay is largely uncharacterized with adolescent inocula, and all storage/freeze–thaw/inoculum-preservation evidence is adult-derived. Novelty here is partly a consequence of missing methodological validation, and should be stated as such.

---

## 9. Recommended next actions

1. **Get the fermentation pH protocol from the bench team** (§3.2) — highest-value single fact; four reports depend on it.
2. **Revisit the butyrate-primacy rationale** in light of §2.2 — the anchor's butyrate contrast was null, so prespecification must carry that choice.
3. **Audit any existing draft text for the *Fusicatenibacter* butyrate association** (§2.1) — this is the finding least likely to survive a specialist reviewer.
4. **Correct the Gillen 2021 citation** from preprint to published Abbott-affiliated trial (§2.4).
5. **Fix the three version traps** before the reference list is built (§4).
6. **Resolve the eight `⚠ UNVERIFIED` items** (§5).
7. **Confirm betadisper was run at 48 h** (§6).
