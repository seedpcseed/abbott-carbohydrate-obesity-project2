# Phase 0 completion — unit inventory

**Date:** 2026-07-19  
**Baseline git SHA:** `67bd320` (committed and pushed before rebuild)  
**Plan:** [`docs/unit-hierarchy-reanalysis-plan.2026-07-19.md`](../docs/unit-hierarchy-reanalysis-plan.2026-07-19.md)

## Acceptance (Phase 0)

| Criterion | Status |
|---|---|
| One table shows 16 SCFA donors, 32 aliquots, 20 planned donors | Met — see `phase0_unit_inventory_summary.csv` |
| Contaminated outputs listed by path | Met — `phase0_archive_do_not_cite.csv`, `phase0_void_claim_cells.csv` |
| No scientific model recompute yet at Phase 0 freeze | Met for inventory; Phase 1+ begins after this file |

## Inventory summary

| Metric | Value |
|---|---|
| Planned donors | 20 |
| Planned aliquots | 40 |
| SCFA analyzed donors | 16 |
| SCFA analyzed aliquots | 32 |
| Missing planned donors | 4 (`19`, `90`, `252`, `287`) |
| Missing planned aliquots | 8 |
| Missing SCFA design samples (48 h) | 4 (`21a nc 48h`, `123a nc 48h`, `66b r2 48h`, `66b s2 48h`) |
| Microbiome rows | 297 |
| Microbiome donors / aliquots | 16 / 32 |
| Donors in both SCFA and 16S | 16 |

## Files in this directory

- `canonical_experimental_units.csv` — SCFA + 16S rows with `donor_id`, `aliquot_id`, `well_id`
- `canonical_experimental_units_scfa.csv`
- `canonical_experimental_units_16s.csv`
- `phase0_unit_inventory_summary.csv`
- `phase0_planned_aliquots.csv`
- `phase0_scfa_analyzed_aliquots.csv`
- `phase0_scfa_missing_design_samples.csv`
- `phase0_scfa_missing_planned_donors.csv`
- `phase0_void_claim_cells.csv`
- `phase0_archive_do_not_cite.csv`

## Frozen SCFA authority (not voided)

`scfa_metabolomics/results/obesity_group_scfa_*.csv` and forest plot from nested Python analysis (N=16 donors) remain the baseline obesity-group SCFA contrasts.
