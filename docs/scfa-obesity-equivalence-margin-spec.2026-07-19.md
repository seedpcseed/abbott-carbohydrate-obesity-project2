# SCFA obesity-group equivalence margin specification

**Date:** July 19, 2026  
**Status:** No primary equivalence margin is locked; formal equivalence is not currently evaluable.

## Estimand and analysis units

The primary estimand is the obesity-minus-healthy-weight difference in 0-to-48-hour concentration change:

`(Obesity 48 h − Obesity 0 h) − (Healthy-weight 48 h − Healthy-weight 0 h)`

It is estimated separately for RDC and SDC for acetate, propionate, and butyrate, giving six primary contrasts. No-added-carbohydrate and the two secondary analytes are sensitivity analyses.

The submission metadata contain 40 A/B aliquot labels from 20 numeric subjects (each person contributes A and B sample splits). R1/R2 and S1/S2 identify biological culture wells within each aliquot; no-added-carbohydrate has one well. The analyzed SCFA file contains 16 of those subjects (8 per obesity group) and 32 aliquots.

## Unit and assay evidence

The DFI-HMMF method documents authentic-standard calibration for acetate, propionate, butyrate, and succinate and describes project QC procedures, including internal-standard recovery, retention time, percent CV, method blanks, and pooled QC. It does not provide this project’s analyte-specific CVs, limits of quantification, total error, or a statement confirming the concentration unit in the exported CSV.

The analysis therefore retains the existing µM label but marks it as pending facility confirmation. Assay precision is treated as a feasibility check, not as a substitute for a biological smallest effect size of interest.

## External evidence review

According to PubMed:

- He et al. described a PFBBr GC-MS SCFA method with detection limits of 0.244–0.977 µM and recoveries of 55.7%–97.9%, but this method and mouse-feces matrix do not establish a biologically negligible obesity-group difference for the present culture system ([DOI](https://doi.org/10.1016/j.jchromb.2018.06.028)).
- Lenzi et al. reported 10%–23% intra/inter-day precision for a different PFBBr method in saliva. That analytical precision is not a biological equivalence threshold for ex vivo fecal fermentation ([DOI](https://doi.org/10.1016/j.jchromb.2023.123826)).
- Comparable in vitro studies in pediatric obesity, normal-weight versus overweight/obese donors, and sorghum fermentation reported SCFA responses or nonsignificant group differences but did not supply reusable equivalence margins for this protocol ([DOI](https://doi.org/10.1128/mBio.00914-20); [DOI](https://doi.org/10.3390/nu13062052); [DOI](https://doi.org/10.3390/nu11020217)).
- Methodological literature emphasizes that the margin determines the conclusion and should be justified using clinical or biological judgment plus statistical reasoning, not selected from the observed treatment difference ([DOI](https://doi.org/10.1136/bmjopen-2024-089587)).

No external source matched the substrate, culture conditions, concentration scale, population, and estimand closely enough to define an analyte-specific negligible difference.

## Locked decision for this analysis version

Version `v0.1-unlocked` records blank lower and upper margins in `scfa_metabolomics/equivalence_margins.csv`.

Consequences:

1. Report model estimates with 95% confidence intervals.
2. Report 90% confidence intervals because they show the interval that a future two-one-sided-test procedure would compare with fixed margins.
3. Report the smallest symmetric bound that would contain each observed 90% interval only as a post hoc precision frontier.
4. Do not report a TOST pass/fail result, non-inferiority, equivalence, preserved capacity, or unaffected capacity.
5. Retain “did not differ significantly” and state that equivalence was not evaluable because an externally justified margin was unavailable.

## Requirements to lock version 1.0

Before adding nonblank margins:

- obtain the facility’s project-level CV, LOQ, dilution, and exported-unit documentation;
- obtain a study-team judgment for the largest negligible 0-to-48-hour difference for each primary SCFA;
- document any directly comparable ex vivo source used to anchor that judgment;
- approve and date the margins without consulting the observed obesity-group contrast estimates;
- rerun `scfa_metabolomics/obesity_equivalence_analysis.py`.

If all six primary 90% confidence intervals then lie inside their prespecified bounds and both one-sided tests have p<0.05, a global equivalence claim may be considered. Otherwise, report analyte/condition-specific results and do not make a global preserved-capacity claim.
