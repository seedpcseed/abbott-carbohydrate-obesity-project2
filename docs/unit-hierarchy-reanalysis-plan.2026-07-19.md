# Unit-hierarchy clean reanalysis plan

**Project:** Abbott–Lurie ZC07, Project 2  
**Date:** July 19, 2026  
**Status:** Executing / Phase 3 refresh underway — baseline `67bd320`; Phase 0 inventory complete; donor-aware integrated HTML re-knit completed 2026-07-19 (`render/integrated-scfa-microbiome.html`).  
**Related:** [`docs/gut-microbes-claim-audit.2026-07-19.md`](gut-microbes-claim-audit.2026-07-19.md), [`docs/gut-microbes-manuscript-outline.2026-07-19.md`](gut-microbes-manuscript-outline.2026-07-19.md), [`docs/scfa-obesity-equivalence-margin-spec.2026-07-19.md`](scfa-obesity-equivalence-margin-spec.2026-07-19.md)

## 1. Decision

The biological unit hierarchy was wrong in most microbiome-facing analyses and in an earlier SCFA “A/B = independent subject” patch. Incremental strata/random-effect edits on top of the current integrated HTML are **not** acceptable.

**Action:** rebuild analysis identity from a single authoritative metadata table, re-specify models once, re-run contaminated analyses end-to-end, then refresh claim-locked documents from new outputs only.

**Do not start Phase 2+ until Phase 0 acceptance criteria are met.**

## 2. Locked biological units

Confirmed by the study team (2026-07-19):

| Level | Label example | Meaning | Independent for group inference? |
|---|---|---|---|
| Subject / donor | `84` | One person (lean or obese) | **Yes** — primary clustering unit |
| Aliquot | `84A`, `84B` | Sample split of that person’s stool | No — nested within subject |
| Culture well | `r1`/`r2`, `s1`/`s2`, `nc` | Replicate culture within aliquot × carbohydrate | No — nested within aliquot |
| Time | `0h`, `48h` | Observation within well | Repeated measure |

**Counts in the SCFA export (current):**

- 16 independent subjects (8 healthy-weight, 8 obesity)
- 32 A/B aliquots
- 160 aliquot × carbohydrate × well cells (design target)
- 4 missing 48 h sample IDs (`21a nc 48h`, `123a nc 48h`, `66b r2 48h`, `66b s2 48h`)

**Counts in submission metadata:**

- 20 numeric subjects / 40 A/B aliquot labels planned
- 4 subjects (8 aliquots) absent from the SCFA quantification export

**Forbidden interpretations going forward:**

- Treating `84A` and `84B` as two people
- Treating culture wells as independent donors
- Calling aliquot×carbohydrate or well rows “N subjects”
- Patching only obesity CIs while leaving microbiome RE/strata on A/B

## 3. Why a clean rebuild (not a patch)

| Area | Current state (2026-07-19) | Contaminated for manuscript? |
|---|---|---|
| Python obesity equivalence (`scfa_metabolomics/obesity_equivalence_analysis.py`) | Nested subject → aliquot → well; N=16 | **No** — retain as SCFA group-contrast authority |
| Integrated SCFA `lmer` source parsing | Updated to nested IDs | Source OK; **HTML not claim-locked until re-knit after microbiome rebuild** |
| Microbiome `subject_id` in phyloseq | Still A/B label (`84B`) | **Yes** |
| Alpha Wilsoxons | Sample/well rows treated as independent | **Yes** |
| PERMANOVA | `strata = subject_id` (A/B) | **Yes** |
| ANCOM-BC2 / MaAsLin2 | Random effect = A/B | **Yes** |
| Taxon–SCFA models | Same RE problem; joins on A/B are fine as keys | **Yes** (inference) |
| Responder overlap / DA / Fusicatenibacter tests | Built on A/B or sample independence | **Yes** |
| Claim audit community/responder numbers | From contaminated models | **Yes** — rewrite after rerun |
| Legacy SCFA Rmds (`*improved*`, `scfa-project2-analysis.Rmd`) | Collapse-average numeric IDs | Archive only — do not cite |

Patching `strata` or `rand_formula` strings inside the existing 6k-line Rmd without freezing metadata and design rules will reintroduce the same error in joins, responder denominators, and claim tables.

## 4. Reusable vs rebuild

### Reusable without recompute

- Raw SCFA quantification CSV and submission-sheet sample IDs
- 16S BIOM / ASV table, taxonomy, qPCR absolute-abundance inputs
- Customer-label strings for matching (`84B.R1.48H`-style)
- Python obesity contrast CSVs / forest plot (after Phase 0 affirms them as frozen SCFA outputs)
- Facility method docs already filed under `docs/` / `scfa_metabolomics/`

### Must rebuild once

1. Authoritative experimental-unit metadata (Phase 0)
2. Microbiome `sample_data` columns and all RE/strata definitions (Phase 1–2)
3. Alpha, beta, DA, taxon–SCFA, responder inference (Phase 2)
4. Integrated HTML/PDF render (Phase 3)
5. Claim audit + outline numerical cells that depend on microbiome/responders (Phase 3)

### Explicitly out of scope for this rebuild

- Locking TOST equivalence margins (still unlocked; separate margin spec)
- Facility confirmation of µM / CV / LOQ (blocking for claim language, not for unit hierarchy)
- Mechanistic *Fusicatenibacter* validation experiments
- Editing legacy SCFA Rmds other than marking them historical if touched

## 5. Target design specification

### 5.1 Canonical metadata columns

Create **one** table (CSV + code object) used by SCFA and microbiome joins:

| Column | Definition |
|---|---|
| `sample_id` | Facility / export sample string (e.g. `84b r2 48h` or microbiome SampleID) |
| `donor_id` | Numeric person (`84`) — **primary independent unit** |
| `aliquot_id` | `84A` / `84B` — nest under donor |
| `carbohydrate` | `No Added Carb` / `RDC` / `SDC` (harmonized labels) |
| `well_repeat` | `1`, `2`, or `1` for NC |
| `well_id` | `donor_id × aliquot_id × carbohydrate × well_repeat` |
| `timepoint` | `0H` / `48H` |
| `group` | healthy-weight / obesity (single mapping table) |
| `assay` | `scfa` / `16S` / `qpcr` as needed |

Demote current microbiome `subject_id` to `aliquot_id`. Never use `aliquot_id` as the top-level random effect for person-level questions.

### 5.2 Model templates (prespecify before coding)

**SCFA concentration (primary structure already aligned in Python / Rmd source):**

```text
concentration ~ group * carbohydrate * timepoint
  + (1 | donor_id)
  + (1 | aliquot_id)
  + (1 | well_id)
```

**Obesity-group estimand (unchanged; keep Python as authority):**

```text
(Obesity 48h − Obesity 0h) − (Healthy-weight 48h − Healthy-weight 0h)
```

Separately for acetate / propionate / butyrate × RDC / SDC. Equivalence remains not evaluable until margins lock.

**Microbiome alpha (replace Wilsoxons for claim-facing tests):**

- Prefer donor-aware mixed models on diversity metrics, **or**
- Analyze one aggregated value per `donor_id × carbohydrate × timepoint` (mean across aliquots/wells) with methods that match that reduced table
- Do not claim person-level effects from well-level Mann–Whitney tests

**PERMANOVA:**

- Distances computed at the chosen observational grain (document: sample vs aggregated donor-cell)
- `strata` / blocks = `donor_id` when multiple rows per donor remain
- Export **term-level** R² and p (not only whole-model summaries)
- Report multivariate dispersion diagnostics

**Differential abundance / taxon–SCFA (ANCOM-BC2, MaAsLin2):**

```text
random effects: (1 | donor_id)
optional nesting: (1 | aliquot_id) if the method and sample size support it
```

If nesting is unstable, collapse to one row per `donor_id × carbohydrate × timepoint` **before** fitting and document the collapse rule.

**Responders (exploratory only):**

- Primary definition: median split of ΔSCFA within analyte × carbohydrate at **`donor_id`** grain (average aliquots/wells first)
- Secondary sensitivity: aliquot-level labels (for matching existing A/B microbiome rows)
- Downstream composition tests must use `donor_id` strata/RE
- Never label record counts as “N subjects”

## 6. Phased execution

### Phase 0 — Freeze units and inventory (no scientific recompute)

**Deliverables:**

1. This plan accepted (checkbox at end).
2. Written unit dictionary matching Section 2.
3. Inventory CSV of: planned donors, analyzed SCFA donors, missing SCFA sample IDs, microbiome aliquot IDs vs `donor_id`.
4. Explicit “archive / do not cite” list for legacy HTML/PDF and old claim numbers.

**Acceptance criteria:**

- [ ] One table shows 16 SCFA donors, 32 aliquots, 20 planned donors
- [ ] No analysis script is edited yet except scaffolding for the metadata table if needed to generate the inventory
- [ ] Contaminated outputs listed by path

### Phase 1 — Metadata rewrite (code, no claim numbers yet)

**Deliverables:**

1. Canonical metadata builder used by integrated Rmd (and optionally mirrored for Python SCFA).
2. Phyloseq `sample_data` gains `donor_id`, `aliquot_id`, `well_id`; A/B is no longer the primary subject.
3. SCFA–microbiome join keys documented: match on `aliquot_id` (+ time + carb); cluster inference on `donor_id`.
4. Knitr cache strategy: **purge** integrated microbiome-related caches after metadata change.

**Acceptance criteria:**

- [ ] `n_distinct(donor_id)` in microbiome sample_data equals the number of numeric persons represented (≈16 overlapping with SCFA, plus any microbiome-only if present — report both)
- [ ] Spot-check: sample `84B.R1.48H` → `donor_id=84`, `aliquot_id=84B`, well=1, time=48H
- [ ] No RE/strata still named solely on A/B without also documenting donor nesting

### Phase 2 — Full scientific rerun from corrected units

Run in this order; each step writes versioned outputs under a new dated results folder (do not overwrite contamination silently):

1. SCFA nested `lmer` + carbohydrate pairwise change contrasts (refresh Tukey SDC-vs-RDC butyrate)
2. Confirm Python obesity contrasts still match the nested R structure (reuse existing CSV if identical)
3. Alpha diversity — donor-aware
4. PERMANOVA — term-level + strata on `donor_id`
5. ANCOM-BC2 / MaAsLin2 — carbohydrate, group, time as prespecified; RE on `donor_id`
6. Taxon–SCFA associations
7. Responder definitions + overlap counts with correct denominators
8. Responder community / DA / focal-taxon tests (secondary, exploratory)

**Acceptance criteria:**

- [ ] Every claim-facing table states the independent unit (`donor_id` N)
- [ ] No Wilcoxon claim-facing result uses unaggregated wells
- [ ] Responder “N” language never says “subjects” for record counts
- [ ] Diagnostics: singularity/boundary notes recorded for nested RE

### Phase 3 — Reporting refresh

**Deliverables:**

1. Re-knit `integrated/integrated-scfa-microbiome.Rmd` → new HTML (and PDF if used)
2. Update claim audit: replace all community/responder numerical cells; keep Python obesity CIs unless Phase 2 shows a change
3. Update manuscript outline Methods/Results unit language and placeholders
4. Mark legacy renders and pre-rebuild claim text as superseded
5. Update sample canvas if inventory counts change ([`scfa-sample-ab-replicates.canvas.tsx`](../../.cursor/projects/home-patcseed-projects-abbott-carbohydrate-obesity-project2/canvases/scfa-sample-ab-replicates.canvas.tsx) — outside repo; optional)

**Acceptance criteria:**

- [ ] Claim audit executive statement still avoids equivalence language unless margins lock
- [ ] HTML Methods describe subject → aliquot → well
- [ ] Diff of claim audit shows which numbers changed and why

### Phase 4 — Stop conditions / residual blockers

These remain **outside** the unit rebuild but still block some claim language:

- Facility confirmation of concentration unit, dilution, analyte CV/LOQ
- Documented reason for four missing SCFA 48 h observations and four absent planned donors
- Externally justified equivalence margins (if desired later)
- Term-level PERMANOVA and pairwise SCFA contrasts must exist before locking those Specific Results sentences

## 7. File ownership during rebuild

| Path | Role in rebuild |
|---|---|
| `docs/unit-hierarchy-reanalysis-plan.2026-07-19.md` | This plan (do not silently diverge) |
| New: `integrated/metadata/` or `microbiome_analysis/metadata/` canonical unit CSV | Phase 0–1 authority |
| `integrated/integrated-scfa-microbiome.Rmd` | Primary rewrite target for microbiome + re-knit |
| `scfa_metabolomics/obesity_equivalence_analysis.py` | Keep; only touch if join columns must rename |
| `scfa_metabolomics/results/obesity_group_scfa_*.csv` | Frozen SCFA group contrasts until a deliberate re-run |
| `microbiome_analysis/analysis.v4_zymo.Rmd` | Secondary; either port unit fix or mark non-manuscript |
| `microbiome_analysis/ANALYSIS_PLAN_PROJECT2.md` | Update subjectID definition after Phase 1 |
| Legacy `scfa_metabolomics/*improved*.Rmd` | Do not cite; optional banner comment only |
| `docs/gut-microbes-claim-audit.2026-07-19.md` | Phase 3 refresh |
| `docs/gut-microbes-manuscript-outline.2026-07-19.md` | Phase 3 refresh |

## 8. Communication rules during rebuild

- Prefer “did not differ significantly” for obesity SCFA; never “equivalent / preserved / unaffected” without locked margins.
- Prefer “associated with” for community findings after rerun.
- Flag every exploratory median-split responder result as secondary.
- If a result vanishes after correct nesting, that is a finding — do not resurrect A/B-as-subject N to recover it.

## 9. Estimated effort (order of work, not calendar commit)

| Phase | Nature |
|---|---|
| 0 | Documentation + inventory (short) |
| 1 | Metadata plumbing + cache purge (moderate; highest leverage) |
| 2 | Full model reruns (longest; especially DA + responders) |
| 3 | Report + claim doc refresh (moderate) |

## 10. Acceptance gate before any recoding

Sign-off checklist:

- [x] Section 2 units accepted as locked
- [x] Clean-rebuild (not patch) approach accepted
- [x] Phase order accepted (0 → 1 → 2 → 3)
- [x] Python obesity contrasts retained as current SCFA group authority
- [x] Contaminated microbiome/responder claim numbers will be treated as void until Phase 3
- [x] Owner ready to execute Phase 0 inventory next

**Approval line:**  

- Accepted by: user (execute after baseline commit/push) Date: 2026-07-19  
- Baseline SHA: `67bd320`  
- Notes / deviations: Phase 0–2 code and metadata landed; full HTML re-knit in progress; claim-audit community cells stay void until Phase 3 refresh.

---

## Appendix A — Contaminated claim cells to void until Phase 3

Treat as **not claim-locked** pending rerun (exact numbers from pre-rebuild audit):

- Healthy–Obese alpha Wilcoxon p / BH values
- Aggregate PERMANOVA R²=0.171 / 0.0467 whole-model summaries (replace with term-level, donor-aware)
- MaAsLin / ANCOM-BC2 coefficients and q-values that used A/B random effects
- Responder alpha / PERMANOVA / absolute *Fusicatenibacter* Wilcoxon p-values and Ns
- Any “N = 48 / 56 / 63 / 94 subjects” style denominators tied to aliquot or well rows

**Retain for now (SCFA person-level, nested Python):**

- Six primary obesity difference-in-change estimates and CIs (acetate/propionate/butyrate × RDC/SDC), with equivalence still not evaluable

## Appendix B — What “done” looks like

The rebuild is complete when:

1. One metadata definition drives SCFA and microbiome.
2. Every primary manuscript analysis clusters people on `donor_id`.
3. Integrated HTML is regenerated after that change.
4. Claim audit community/responder cells are rewritten from the new run (or explicitly retired).
5. No live workflow still documents `84B` as an independent subject for inference.
