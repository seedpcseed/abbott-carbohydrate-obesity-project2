# Media / negative-control handling (provisional)

**Date:** 2026-07-19  
**Status:** Documented from available repository assets; facility QC confirmation still pending.

## What is in the repo

- Study culture design includes a no-added-carbohydrate arm (`No Added Carb` / `nc`) nested under the same donor → aliquot → well hierarchy as RDC/SDC.
- Matched SCFA–16S–qPCR tables retain no-added-carbohydrate wells when both assays are present (`tables/matched_scfa_16s_qpcr_both.csv`).
- Primary carbohydrate-change contrasts and the focal *Fusicatenibacter* estimand use **SDC − no-carb** (and RDC − no-carb for specificity), so the no-carb arm is the within-donor fermentation baseline rather than a discarded blank.

## What is not yet confirmed

- Dedicated media-only or extraction-negative qPCR Ct / copy-number results are not present as a locked QC table.
- Primer target, standard curve, LOD, and LOQ remain facility-blocking (Gate 1).
- No automatic exclusion rule for near-LOD Ct values is applied beyond requiring `gene_copies_per_ul > 0` for qPCR-scaled modeling.

## Rule used in additional analyses

1. Do not drop no-added-carbohydrate culture wells from the matched table.
2. Use them as the paired baseline for net-change estimands.
3. Do not interpret media blanks or negatives without facility files; keep absolute-language claims provisional.
