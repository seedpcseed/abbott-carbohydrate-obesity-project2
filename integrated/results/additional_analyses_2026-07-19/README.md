# Additional analyses results — 2026-07-19

**Plan:** [`docs/additional-analyses-plan.2026-07-19.md`](../../docs/additional-analyses-plan.2026-07-19.md)  
**Status:** Priority A/B executed. MaAsLin 3 / ALDEx3 ran in a separate R 4.6 / Bioconductor 3.23 library under `integrated/additional_analyses/renv_bioc323/`. Facility qPCR LOD/LOQ/primer documentation remains **blocking** for absolute units language. Integrated HTML re-knit includes claim-retirement wording for sample-level *Fusicatenibacter* panels.

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
| Genus matrices for absolute DA | `tables/genus_relative_for_da.csv`, `genus_counts_matched.csv` |
| MaAsLin 3 + ALDEx3 qPCR-scaled DA | `absolute_da/` (see that folder’s README) |

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

## Remaining blockers

1. Facility confirmation of qPCR primers, efficiency, LOD/LOQ, and the `genome_copies_per_ul*` conversion footnote.
2. Figure audit for independent-donor display (Gate 4 residual).
3. Optional full ANCOM-BC2 coefficient export for formal three-method concordance tables (MaAsLin3/ALDEx3 coefficient extracts are in `absolute_da/`).

## Reproduce

```bash
bash integrated/additional_analyses/run_all.sh
```
