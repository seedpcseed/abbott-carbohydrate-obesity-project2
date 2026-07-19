# Additional analyses results — 2026-07-19

**Plan:** [`docs/additional-analyses-plan.2026-07-19.md`](../../docs/additional-analyses-plan.2026-07-19.md)  
**Status:** Priority A/B core analyses executed in the current R 4.3 / Bioconductor 3.18 and Python environments. MaAsLin 3 / ALDEx3 require a separate R ≥4.6 / Bioconductor 3.23 lockfile and are **not** run yet. Facility qPCR LOD/LOQ/primer documentation remains **blocking** for absolute units language.

## What was run

| Step | Output |
|---|---|
| Exact well-matched SCFA–16S–qPCR table | `tables/matched_scfa_16s_qpcr*.csv`, `data_audit/` |
| Sample / missingness reconciliation | `data_audit/phaseA_inventory.csv`, missing-donor/sample CSVs |
| Carbohydrate-change SCFA contrasts (acetate, propionate, butyrate) | `scfa_contrasts/carbohydrate_change_contrasts.csv` |
| Donor-aware alpha models | `alpha_permanova/alpha_*.csv` |
| Term-level PERMANOVA + betadisper (donor×carb×time aggregation) | `alpha_permanova/permanova_*.csv`, `betadisper_48h.csv` |
| Total bacterial load donor contrasts | `total_load/` |
| Focal donor-level qPCR-scaled *Fusicatenibacter*–butyrate | `focal_fusicatenibacter/` |
| Variance components / ICC sketch | `diagnostics/variance_components_deltas.csv` |

## Primary focal result (exploratory)

SDC-attributable Δlog10(qPCR-scaled *Fusicatenibacter*) vs SDC-attributable Δbutyrate:

- Spearman ρ ≈ **0.48**
- Donor permutation p ≈ **0.088**
- 95% donor-bootstrap CI for ρ ≈ **−0.24 to 0.89**
- N = **14** donors with complete net-change pairs
- OLS slope ≈ **1.39 µM per log10 abundance unit** (descriptive)

Direction is positive but the confidence interval includes zero. Do not claim a significant association. Relative-abundance and total-load specificity checks were nonsignificant after BH correction within the secondary family.

## Key SCFA / community outputs

- SDC versus RDC change contrasts remain nonsignificant for acetate, propionate, and butyrate under Sidak-adjusted pairwise change contrasts (see `scfa_contrasts/`).
- 48-hour donor-aggregated PERMANOVA: group R²≈0.055, p≈0.0027; carbohydrate p≈0.41; interaction p≈1.0.
- betadisper at 48 h: carbohydrate p≈0.88; group p≈0.094.
- Alpha mixed-model Healthy-weight vs obesity at 48 h under no added carb: Shannon and Observed raw group contrasts p<0.01 (see `alpha_group_contrasts.csv`); interpret with multiplicity family and nestings as specified in the plan.

## Blockers / not executed

1. Facility confirmation of qPCR primers, efficiency, LOD/LOQ, and the `genome_copies_per_ul*` conversion footnote.
2. Separate locked environment for **MaAsLin 3** and **ALDEx3** (R ≥4.6, Bioconductor 3.23).
3. Full integrated HTML re-knit incorporating these exports (optional next step; results are already versioned under this folder).
4. Mechanistic / targeted *Fusicatenibacter* qPCR (Priority C).

## Reproduce

```bash
bash integrated/additional_analyses/run_all.sh
```
