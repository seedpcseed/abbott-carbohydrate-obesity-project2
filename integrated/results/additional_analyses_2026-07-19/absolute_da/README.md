# Modern absolute / qPCR-scaled DA results

**Date:** 2026-07-19  
**Environment:** R 4.6.0 + Bioconductor 3.23 project library at `integrated/additional_analyses/renv_bioc323/library`  
**Packages:** maaslin3 1.4.0; ALDEx3 1.2.0

## What ran

| Analysis | Key outputs |
|---|---|
| qPCR scale verify | `maaslin3_qpcr_scale_verify.txt` (scaled == rel × gene_copies) |
| MaAsLin 3 qPCR-scaled abundance | `maaslin3_abundance_results.csv`, `all_results.tsv`, `significant_results.tsv` |
| MaAsLin 3 prevalence | `maaslin3_prevalence_results.csv` |
| ALDEx3 sample.sm grid | `aldex3_results_svar_{0.05,0.25,1.00}.csv`, `aldex3_results_all_svar.csv` |
| ALDEx3 exact lme4 trim | `aldex3_results_svar_0.25_lme4_trimmed.csv` |
| Focal taxon extracts | `fusicatenibacter_maaslin3.csv`, `fusicatenibacter_aldex3_svar0.25.csv` |
| ANCOM-BC2 relative abundances | `ancombc2_relative_results_{wide,long}.csv` |
| Three-method coefficient concordance | `concordance_ancombc2_maaslin3_aldex3.csv`, `concordance_focus_effects.csv`, `concordance_fusicatenibacter.csv` |

## Model notes

- Scaling: genus relative × `gene_copies_per_ul`; MaAsLin3 `normalization=NONE`, `transform=LOG`, `median_comparison_abundance=FALSE`.
- Formula used: `~ carbohydrate * timepoint + obesity_group + (1|donor_id)`.
- Dual random effect `(1|donor_id)+(1|aliquot_id)` triggered a maaslin3 1.4.0 post-fit assembly error on this table; documented and skipped.
- ALDEx3: `sample.sm` with `s.mu=log2(gene_copies_per_ul)` and provisional `s.var ∈ {0.05,0.25,1.0}` because facility qPCR technical variance is unavailable; engine `blmm` for full table, `lme4` on a trimmed feature set.
- Counts for ALDEx3 ≈ `round(genus_relative × asv_read_depth)` with >75% zero genera amalgamated to `other`.
- ANCOM-BC2 (relative): render-env ANCOMBC 2.4; `fix_formula = carb_time + obesity_group` with `carb_time = carbohydrate_timepoint` because this build does not parse `*`/`:` interaction strings; `(1|donor_id)` random intercept when available.
- Concordance compares coefficient **directions** on loosely aligned effect labels. Interaction estimands are not identical across packages (MaAsLin/ALDEx interaction terms vs ANCOM carb_time levels vs reference), so discordant signs are diagnostic rather than automatic contradictions.

## Interpretive lock

- These are exploratory genus-wide screens, not a second claim-facing *Fusicatenibacter*–butyrate test.
- *Fusicatenibacter* carbohydrateSDC:timepoint48H is directionally positive in MaAsLin3 and ALDEx3 (s.var=0.25) but not significant after multiplicity (MaAsLin joint q≈0.13; ALDEx adj p≈0.58). The ANCOM-BC2 relative `carb_timeSDC_48H` contrast vs reference is near zero/negative (q≈0.70), consistent with a compositional/relative signal that does not match the qPCR-scaled abundance interaction.
- Across focus effects with all three coefficients available, 85/182 rows are direction-concordant; for SDC×48H specifically 16/45 are three-way concordant.
- Facility LOD/LOQ/primer docs remain blocking for absolute units language.

## Reproduce

```bash
export R_LIBS="$PWD/integrated/additional_analyses/renv_bioc323/library"
python3 integrated/additional_analyses/python/prepare_genus_matrices.py
Rscript integrated/additional_analyses/R/absolute_da_maaslin3.R
Rscript integrated/additional_analyses/R/absolute_da_aldex3.R
# ANCOM-BC2 uses the render env (BioC 3.18 / ANCOMBC 2.4):
env -u R_LIBS -u R_LIBS_USER conda run -p /tmp/abbott-render-r \
  Rscript integrated/additional_analyses/R/absolute_da_ancombc2_concordance.R
```
