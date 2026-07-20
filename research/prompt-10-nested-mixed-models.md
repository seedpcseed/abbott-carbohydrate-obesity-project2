# Prompt 10 — Statistical citations for nested ex vivo fermentation designs (LMM / ICC / PERMANOVA / diagnostics)

**Manuscript target:** Methods 2.8; Results 3.5–3.6; Discussion (strengths/limitations) · **Citation job:** Supply verified methodological and software citations for donor-nested linear mixed models, ICC-based repeatability, PERMANOVA/betadisper on Bray–Curtis distances, and resampling-based sensitivity analysis · **Date:** 2026-07-19

**Verification status:** Every DOI/PMID below was retrieved from PubMed E-utilities or the Crossref REST API during this run. Items marked `⚠ UNVERIFIED` were not confirmed and must be checked manually before submission. (PubMed-sourced records are attributed to PubMed; DOI links are given for each.)

---

## Prioritized citation table

### A. Linear mixed models, nested random effects, estimated marginal means

| # | Title | Year | Journal / Publisher | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 1 | Fitting Linear Mixed-Effects Models Using lme4 | 2015 | J Stat Softw 67(1) | [10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01) | Canonical software citation for the LMM fits; required whenever `lmer()` is used |
| 2 | lmerTest Package: Tests in Linear Mixed Effects Models | 2017 | J Stat Softw 82(13) | [10.18637/jss.v082.i13](https://doi.org/10.18637/jss.v082.i13) | Satterthwaite/Kenward–Roger denominator df and Type III F-tests for fixed effects |
| 3 | emmeans: Estimated Marginal Means, aka Least-Squares Means | 2017– | CRAN Contributed Packages (Lenth & Piaskowski) | [10.32614/CRAN.package.emmeans](https://doi.org/10.32614/CRAN.package.emmeans) | Software citation for EMMs and Tukey-adjusted pairwise contrasts |
| 4 | Mixed Effects Models and Extensions in Ecology with R — esp. Ch. 5, "Mixed Effects Modelling for Nested Data" | 2009 | Springer (Zuur, Ieno, Walker, Saveliev, Smith) | Book [10.1007/978-0-387-87458-6](https://doi.org/10.1007/978-0-387-87458-6); Ch.5 [10.1007/978-0-387-87458-6_5](https://doi.org/10.1007/978-0-387-87458-6_5) · ISBN 978-0-387-87457-9 | **Best single textbook cite for the nested random-effects structure** (donor / aliquot-within-donor / well-within-aliquot) |
| 5 | Mixed-Effects Models in S and S-PLUS | 2000 | Springer (Pinheiro & Bates) | [10.1007/978-1-4419-0318-1](https://doi.org/10.1007/978-1-4419-0318-1) · ISBN 978-1-4419-0317-4 | Foundational treatment of nested/crossed random effects and repeated measures |
| 6 | A brief introduction to mixed effects modelling and multi-model inference in ecology | 2018 | PeerJ 6:e4794 | PMID 29844961 · [10.7717/peerj.4794](https://doi.org/10.7717/peerj.4794) | Accessible best-practice guide; explicitly warns against fitting random effects with too few levels and against treating pseudoreplicates as independent |
| 7 | Data Analysis Using Regression and Multilevel/Hierarchical Models | 2006 | Cambridge Univ. Press (Gelman & Hill) | [10.1017/CBO9780511790942](https://doi.org/10.1017/CBO9780511790942) · ISBN 978-0-521-68689-1 | Standard reference for partial pooling and variance decomposition across hierarchical levels |
| 8 | Variance Components | 1992 | Wiley (Searle, Casella, McCulloch) | [10.1002/9780470316856](https://doi.org/10.1002/9780470316856) · ISBN 978-0-471-62162-1 | Formal statistical authority for REML variance-component estimation underlying the nested decomposition |
| 9 | Evaluating significance in linear mixed-effects models in R | 2017 (online 2016) | Behav Res Methods 49:1494–1502 | [10.3758/s13428-016-0809-y](https://doi.org/10.3758/s13428-016-0809-y) | Justifies Satterthwaite/KR df over naive Wald z for small numbers of clusters (n=16 donors) |
| 10 | A Kenward-Roger Approximation and Parametric Bootstrap Methods for Tests in Linear Mixed Models (pbkrtest) | 2014 | J Stat Softw 59(9) | [10.18637/jss.v059.i09](https://doi.org/10.18637/jss.v059.i09) | Cite if KR df or parametric-bootstrap LRTs are reported |
| 11 | Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing | 1995 | J R Stat Soc B 57:289–300 | [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) | If FDR (rather than Tukey) is used across the three SCFA families |

### B. ICC / repeatability of donor-level responses

| # | Title | Year | Journal | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 12 | A Guideline of Selecting and Reporting Intraclass Correlation Coefficients for Reliability Research | 2016 | J Chiropr Med 15:155–163 | PMID 27330520 · [10.1016/j.jcm.2016.02.012](https://doi.org/10.1016/j.jcm.2016.02.012) | **Primary ICC cite.** Specifies that ICC form (model/type/definition) and software must be reported; supplies the poor/moderate/good/excellent <0.5 / 0.5–0.75 / 0.75–0.9 / >0.9 bands |
| 13 | Repeatability for Gaussian and non-Gaussian data: a practical guide for biologists | 2010 | Biol Rev 85:935–956 | PMID 20569253 · [10.1111/j.1469-185X.2010.00141.x](https://doi.org/10.1111/j.1469-185X.2010.00141.x) | **Primary repeatability cite.** Defines ICC as the between-subject share of phenotypic variance and endorses LMM-based estimation with parametric-bootstrap CIs — exactly the donor-repeatability framing needed |
| 14 | rptR: Repeatability Estimation for Gaussian and Non-Gaussian Data | 2016– | CRAN Contributed Packages (Stoffel, Nakagawa, Schielzeth) | [10.32614/CRAN.package.rptR](https://doi.org/10.32614/CRAN.package.rptR) | Software cite if rptR is used for bootstrap ICC confidence intervals |
| 15 | Within- and between-child variation in repeated urinary pesticide metabolite measurements over a 1-year period | 2014 | Environ Health Perspect 122:201–206 | PMID 24325925 · [10.1289/ehp.1306737](https://doi.org/10.1289/ehp.1306737) | *Backup only.* Worked precedent for LMM-derived ICCs partitioning within- vs between-person biomarker variance; different matrix (urinary pesticides), cite for method not for topic |

### C. PERMANOVA, multivariate dispersion, and Bray–Curtis

| # | Title | Year | Journal / Publisher | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 16 | A new method for non-parametric multivariate analysis of variance | 2001 | Austral Ecol 26:32–46 | [10.1111/j.1442-9993.2001.01070.pp.x](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x) | **Original PERMANOVA method cite** |
| 17 | Fitting multivariate models to community data: a comment on distance-based redundancy analysis | 2001 | Ecology 82:290–297 | [10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2](https://doi.org/10.1890/0012-9658(2001)082%5B0290:FMMTCD%5D2.0.CO;2) | McArdle & Anderson — the partitioning of a dissimilarity matrix that gives **term-specific sums of squares and R²**; the citation to use when reporting per-term R² rather than an omnibus result |
| 18 | Distance-based tests for homogeneity of multivariate dispersions | 2006 | Biometrics 62:245–253 | PMID 16542252 · [10.1111/j.1541-0420.2005.00440.x](https://doi.org/10.1111/j.1541-0420.2005.00440.x) | **Original betadisper method cite** (PERMDISP); multivariate Levene's test on principal-coordinate distances to centroid/spatial median |
| 19 | PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? | 2013 | Ecol Monogr 83:557–574 | [10.1890/12-2010.1](https://doi.org/10.1890/12-2010.1) | **The location-vs-dispersion caveat.** Shows a significant PERMANOVA can reflect dispersion heterogeneity rather than centroid shift, especially in unbalanced designs — mandatory alongside any betadisper result |
| 20 | Distance-based multivariate analyses confound location and dispersion effects | 2012 (online 2011) | Methods Ecol Evol 3:89–101 | [10.1111/j.2041-210X.2011.00127.x](https://doi.org/10.1111/j.2041-210X.2011.00127.x) | Warton, Wright & Wang — independent statement of the same confounding, framed via the mean–variance relationship of abundance data |
| 21 | Permutational Multivariate Analysis of Variance (PERMANOVA) | 2017 | Wiley StatsRef: Statistics Reference Online | [10.1002/9781118445112.stat07841](https://doi.org/10.1002/9781118445112.stat07841) | Anderson's own modern review; covers restricted permutation for nested/blocked designs — the correct cite for permuting **within donor** |
| 22 | Measures of precision for dissimilarity-based multivariate analysis of ecological communities | 2015 (online 2014) | Ecol Lett 18:66–73 | [10.1111/ele.12385](https://doi.org/10.1111/ele.12385) | Anderson & Santana-Garcon; effect-size/precision reporting for dissimilarity-based analyses beyond a bare p-value |
| 23 | vegan: Community Ecology Package | 2001– | CRAN Contributed Packages (Oksanen, Simpson, Blanchet, Kindt, Legendre *et al.*) | [10.32614/CRAN.package.vegan](https://doi.org/10.32614/CRAN.package.vegan) | Software cite for `adonis2()`, `betadisper()`, `vegdist()`; state the package version |
| 24 | An Ordination of the Upland Forest Communities of Southern Wisconsin | 1957 | Ecol Monogr 27:325–349 | [10.2307/1942268](https://doi.org/10.2307/1942268) | Bray & Curtis — origin of the Bray–Curtis dissimilarity |
| 25 | Numerical Ecology (3rd English ed.) | 2012 | Elsevier, Developments in Environmental Modelling 24 (Legendre & Legendre) | [10.1016/C2010-0-66470-4](https://doi.org/10.1016/C2010-0-66470-4) · ISBN 978-0-444-53868-0 | Textbook authority for dissimilarity-coefficient choice and permutation testing in community ecology |

### D. Nested replicates, experimental unit, and design precedents

| # | Title | Year | Journal | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 26 | Pseudoreplication and the Design of Ecological Field Experiments | 1984 | Ecol Monogr 54:187–211 | [10.2307/1942661](https://doi.org/10.2307/1942661) | Hurlbert — the canonical statement that inference must be at the level of the true experimental unit (here: the donor) |
| 27 | The problem of pseudoreplication in neuroscientific studies: is it affecting your analysis? | 2010 | BMC Neurosci 11:5 | [10.1186/1471-2202-11-5](https://doi.org/10.1186/1471-2202-11-5) | Lazic — **explicitly endorses averaging technical/sub-samples to the biological unit** as a valid alternative to a hierarchical model; the cleanest precedent for collapsing wells before donor-level microbiome inference |
| 28 | Best practices for analysing microbiomes | 2018 | Nat Rev Microbiol 16:410–422 | PMID 29795328 · [10.1038/s41579-018-0029-9](https://doi.org/10.1038/s41579-018-0029-9) | Knight *et al.* — general design/reproducibility authority for microbiome studies |
| 29 | The Power of Microbiome Studies: Some Considerations on Which Alpha and Beta Metrics to Use and How to Report Results | 2022 | Front Microbiol 12:796025 | [10.3389/fmicb.2021.796025](https://doi.org/10.3389/fmicb.2021.796025) | Kers & Saccenti — reporting standards for beta-diversity metrics, effect size, and power; supports reporting R² alongside p |
| 30 | Analysing microbiome intervention design studies: Comparison of alternative multivariate statistical methods | 2021 | PLoS ONE 16:e0259973 | [10.1371/journal.pone.0259973](https://doi.org/10.1371/journal.pone.0259973) | Khomich *et al.* — benchmark showing community-level tests (incl. PERMANOVA) agree on effect size/significance while taxon-ranking methods diverge; supports treating community-level results as the robust layer |
| 31 | Specific substrate-driven changes in human faecal microbiota composition contrast with functional redundancy in short-chain fatty acid production | 2018 | ISME J 12:610–622 | PMID 29192904 · [10.1038/ismej.2017.196](https://doi.org/10.1038/ismej.2017.196) | **Design precedent only:** in vitro batch fermentation with donor as the biological unit and carbohydrate substrate as the within-donor factor. See "Conflicts / nuance" — do **not** cite for functional-redundancy conclusions |

### E. Diagnostics, robustness, and resampling sensitivity analysis

| # | Title | Year | Journal / Publisher | DOI/PMID | Why cite |
|---|---|---|---|---|---|
| 32 | Robustness of linear mixed-effects models to violations of distributional assumptions | 2020 | Methods Ecol Evol 11:1141–1152 | [10.1111/2041-210X.13434](https://doi.org/10.1111/2041-210X.13434) | Schielzeth *et al.* — **the key diagnostics cite.** Simulation evidence that LMM fixed-effect estimates, SEs and variance components are largely unbiased under non-normal and heteroscedastic residuals; justifies retaining the LMM and reporting resampling as sensitivity rather than switching models |
| 33 | An Introduction to the Bootstrap | 1993/1994 | Chapman & Hall/CRC (Efron & Tibshirani) | [10.1201/9780429246593](https://doi.org/10.1201/9780429246593) · ISBN 978-0-412-04231-7 | Standard bootstrap authority for resampling-based CIs |
| 34 | Randomization, Bootstrap and Monte Carlo Methods in Biology | 2018 (3rd/4th ed.) | Chapman & Hall/CRC (Manly) | [10.1201/9781315273075](https://doi.org/10.1201/9781315273075) · ISBN 978-1-315-27307-5 | Permutation-test authority for biological data; pairs with #21 for restricted (within-donor) permutation |

> Also usable for resampling: **#13 (Nakagawa & Schielzeth 2010)**, which specifically recommends parametric bootstrapping and randomisation tests as the most accurate ways to obtain CIs and p-values for ICC/repeatability estimates.

---

## Prior-claim bullets (mapped to citation numbers)

- **The donor is the experimental unit; aliquots and wells are sub-samples.** Analysing sub-samples as independent replicates inflates the effective sample size and Type I error; inference must be returned to the donor level either by nested random effects or by averaging to the donor. [26, 27, 6]
- **Nested random effects are the standard device for donor / aliquot / well hierarchies.** Random intercepts for donor, aliquot-within-donor, and well-within-aliquot absorb the correlation among repeated timepoints sharing a well. [4, 5, 8, 7]
- **REML variance components partition total variance across hierarchical levels**, which is what licenses a variance-decomposition narrative. [8, 4, 1]
- **Denominator degrees of freedom matter when clusters are few.** With 16 donors, Satterthwaite or Kenward–Roger approximations are preferred over asymptotic Wald tests. [2, 9, 10]
- **Estimated marginal means with Tukey adjustment are the standard route to condition contrasts from an LMM**, and multiplicity must be controlled within a defined family. [3, 11]
- **ICC = between-donor variance ÷ total variance, and is the standard index of interpersonal repeatability.** Reporting requires stating the ICC form and the software, and interpreting against published bands. [12, 13]
- **LMM-based ICC estimation is preferred over ANOVA/correlation approaches** because covariates and unbalanced designs are handled naturally; CIs are best obtained by parametric bootstrap. [13, 14]
- **PERMANOVA tests differences in multivariate centroids using permutation on a dissimilarity matrix**, with sums of squares partitioned directly from the distance matrix. [16, 17, 25]
- **Term-specific R² is the informative quantity, not an omnibus p-value.** McArdle & Anderson's partitioning yields per-term SS/R²; marginal (`by = "margin"`) rather than sequential (`by = "terms"`) testing avoids order-dependence of terms. [17, 23, 29]
- **A significant PERMANOVA does not by itself establish a location (centroid) difference** — heterogeneity of multivariate dispersion can produce the same rejection, particularly in unbalanced designs. betadisper must be reported alongside. [19, 20, 18]
- **Restricted permutation within donor is the correct null for a nested/blocked design.** Free permutation across donors would test the wrong hypothesis. [21, 16]
- **Effect sizes and precision measures should accompany dissimilarity-based tests**, not p-values alone. [22, 29]
- **Averaging technical replicates to the biological unit before analysis is a recognised, defensible alternative to a full hierarchical model.** [27, 26]
- **LMMs are robust to moderate departures from normality and homoscedasticity**; where residual diagnostics are imperfect, permutation/bootstrap sensitivity analyses are the appropriate confirmation, not a change of estimator. [32, 33, 34, 13]

---

## Conflicts / nuance

1. **Reichardt et al. 2018 (#31) is a double-edged citation.** Its title and conclusions assert *functional redundancy* in SCFA production across donors despite compositional differences. That is precisely the shape of claim the project guardrails forbid ("preserved fermentation capacity", obesity-group equivalence). **Cite #31 only as a design precedent** — donor as biological unit, substrate as within-donor factor, n=3 donors — and never as evidence bearing on group equivalence or capacity preservation in this study.
2. **betadisper and PERMANOVA are not a clean two-step decision rule.** Anderson & Walsh (#19) show that a *non-significant* betadisper does not license an unqualified location interpretation, and a *significant* betadisper does not automatically invalidate PERMANOVA. Both are informative; neither is a gate. Report both and interpret jointly.
3. **Warton et al. (#20) partly disagree with Anderson's framing.** They attribute the confounding to the mean–variance relationship of abundance data and argue for model-based multivariate alternatives rather than dissimilarity-based tests. Citing #19 and #20 together honestly signals that the caveat is real and that the field is not unanimous on the remedy.
4. **Sequential vs marginal SS in `adonis2()`.** vegan's historical default was sequential (`by = "terms"`), which makes each term's R² depend on model order. This is a documented and actively debated design choice in the vegan project. State explicitly which `by=` setting was used; #17 is the theoretical backing for the partitioning either way. `⚠` The vegan GitHub discussions on this are developer commentary, not citable literature — do not cite them.
5. **ICC form must be specified.** Koo & Li (#12) enumerate ten ICC forms with different assumptions. An LMM-derived variance-ratio ICC (à la #13) is not identical to Shrout–Fleiss ICC(2,1) as computed by reliability packages. Name the estimator and the variance components in the numerator/denominator.
6. **Koo & Li's benchmark bands are from clinical reliability research**, not microbial ecology. They are widely used but are conventions, not thresholds with a biological basis. Present them as descriptive anchors.
7. **Random effects with few levels.** Harrison et al. (#6) caution against random effects with fewer than ~5 levels. With 16 donors the donor term is fine; an aliquot-within-donor term with 2 levels per donor is at the boundary and may yield near-zero or singular variance estimates. Report singular fits honestly rather than dropping terms silently.
8. **Schielzeth et al. (#32) reduce, but do not eliminate, the need for diagnostics.** Their result is that *fixed-effect* inference is robust; *random-effect* variance estimates and their CIs are more sensitive. Since the manuscript's variance decomposition and ICCs are random-effect quantities, the robustness claim should not be over-extended to them — this is exactly where the bootstrap sensitivity analysis earns its place.
9. **Shrout & Fleiss 1979** (the other standard ICC reference) was **not verified in this run**. `⚠ UNVERIFIED — needs manual check` if it is added to the reference list.
10. **Stoffel, Nakagawa & Schielzeth 2017, *Methods Ecol Evol* 8:1639–1644** (the rptR journal article, as distinct from the CRAN package record #14) was **not verified in this run**. `⚠ UNVERIFIED — needs manual check`.

---

## What we must not overclaim (given claim lock)

- Nothing in this citation set supports **SDC superiority over RDC** for any SCFA. These are inference-machinery citations only; a Tukey-adjusted contrast citation (#3, #11) licenses the *procedure*, never the *direction* of a result.
- A **non-significant group × carbohydrate term does not demonstrate obesity-group equivalence or "preserved fermentation capacity."** Absence of evidence is not equivalence; no citation here supports an equivalence claim, and none of #1–#34 is an equivalence-testing (TOST) reference. Do not let #31's "functional redundancy" language leak into this manuscript.
- A **significant PERMANOVA term does not establish carbohydrate-driven 48-hour community restructuring.** Given #19 and #20, it establishes a difference in multivariate *location and/or dispersion*, and the accompanying R² must be reported so readers can see the magnitude.
- **ICC quantifies repeatability of donor rank/level across conditions. It is not a responder phenotype, not a predictive model, and not evidence that the microbiome predicts dietary response.** High ICC means donors are consistently ordered; it says nothing about whether that ordering is predictable from baseline composition.
- **Term-specific R² is a descriptive variance share, not an effect size with clinical meaning.** Do not translate R² into host-level, glycemic, or obesity-treatment relevance.
- **Bootstrap/permutation agreement is a robustness statement about the analysis, not corroboration of a biological claim.** Phrase it as "conclusions were unchanged under resampling", never as "confirmed".
- These are **statistical** citations. None of them is evidence about *Fusicatenibacter*, butyrate producers, or carbohydrate digestion kinetics.

---

## Downstream synthesis blocks

### Block 1 — Methods 2.8: nested LMM structure

**SECTION:** Methods 2.8 (Statistical analysis — longitudinal SCFA models)
**PRIOR CLAIM:** Where paired aliquots and multiple culture wells are nested within donor and wells link repeated timepoints, the donor is the experimental unit and the hierarchy must be modelled with nested random effects fitted by REML.
**CITATIONS (lead):** #4 (Zuur et al. 2009, Ch. 5) · #5 (Pinheiro & Bates 2000) · #6 (Harrison et al. 2018) — **(backup):** #8 (Searle et al. 1992) · #26 (Hurlbert 1984)
**CONFLICTS / NUANCE:** An aliquot-within-donor term with only two levels per donor sits at the lower bound of what #6 recommends and may return a singular or near-zero variance estimate; this should be reported rather than resolved by silent term removal.
**MANUSCRIPT-SAFE SENTENCE:** "Longitudinal SCFA concentrations were analysed with linear mixed-effects models incorporating nested random intercepts for donor, aliquot within donor, and culture well within aliquot × condition, reflecting the hierarchical structure of the fermentation design and treating the donor as the independent biological unit."
**CLAIM-LOCK CHECK:** Supported (methodological description only).

### Block 2 — Methods 2.8: software, degrees of freedom, and contrasts

**SECTION:** Methods 2.8 (Software and inference)
**PRIOR CLAIM:** Models were fitted with lme4, Type III tests obtained via lmerTest with Satterthwaite denominator degrees of freedom, and condition comparisons drawn as Tukey-adjusted estimated marginal means.
**CITATIONS (lead):** #1 (lme4) · #2 (lmerTest) · #3 (emmeans) — **(backup):** #9 (Luke 2017) · #10 (pbkrtest)
**CONFLICTS / NUANCE:** With 16 donors, asymptotic Wald tests are anticonservative; #9 supports the Satterthwaite/KR choice explicitly. State the multiple-comparison family (within SCFA, across conditions) so the Tukey adjustment is interpretable.
**MANUSCRIPT-SAFE SENTENCE:** "Models were fitted by restricted maximum likelihood using lme4; Type III F-tests with Satterthwaite denominator degrees of freedom were obtained with lmerTest; and estimated marginal means with Tukey-adjusted pairwise contrasts were computed with emmeans."
**CLAIM-LOCK CHECK:** Supported.

### Block 3 — Results 3.5: ICC / interpersonal repeatability

**SECTION:** Results 3.5 (Donor repeatability of SCFA responses)
**PRIOR CLAIM:** The intraclass correlation coefficient — the between-donor share of total variance — is the standard index for how repeatable a donor's biomarker response is across conditions, and is best estimated from the mixed model with bootstrap confidence intervals.
**CITATIONS (lead):** #12 (Koo & Li 2016) · #13 (Nakagawa & Schielzeth 2010) — **(backup):** #14 (rptR) · #15 (Attfield et al. 2014)
**CONFLICTS / NUANCE:** #12 requires that the ICC form, model, type and definition be stated; an LMM variance-ratio ICC is not interchangeable with Shrout–Fleiss forms. #12's poor/moderate/good/excellent bands originate in clinical reliability research and are conventions, not biologically calibrated cut-offs.
**MANUSCRIPT-SAFE SENTENCE:** "Interpersonal repeatability was summarised as the intraclass correlation coefficient, computed from the mixed-model variance components as the proportion of total variance attributable to between-donor differences, with confidence intervals obtained by parametric bootstrap; the ICC form and estimating software are reported in full."
**CLAIM-LOCK CHECK:** Supported as a descriptive variance statistic. **Unsupported if used beyond evidence:** any framing of ICC as identifying responder phenotypes or as evidence that donor microbiome composition predicts dietary response.

### Block 4 — Methods 2.8 / Results 3.6: PERMANOVA on Bray–Curtis

**SECTION:** Methods 2.8 and Results 3.6 (Community composition)
**PRIOR CLAIM:** Bray–Curtis dissimilarities were analysed by PERMANOVA with permutation restricted within donor, and term-specific R² is reported rather than an omnibus result alone.
**CITATIONS (lead):** #16 (Anderson 2001) · #17 (McArdle & Anderson 2001) · #21 (Anderson 2017, restricted permutation) — **(backup):** #23 (vegan) · #24 (Bray & Curtis 1957) · #25 (Legendre & Legendre 2012)
**CONFLICTS / NUANCE:** `adonis2()` sequential testing makes term R² order-dependent; the `by=` setting and package version must be stated. Term-specific R² should be reported for every term, not only significant ones (#29).
**MANUSCRIPT-SAFE SENTENCE:** "Community-level differences in Bray–Curtis dissimilarity were tested by permutational multivariate analysis of variance (adonis2, vegan v.X.X-X) with permutations restricted within donor to respect the paired design; term-specific R² values are reported for all model terms alongside permutation p-values."
**CLAIM-LOCK CHECK:** Supported. **Unsupported if used beyond evidence:** describing any significant term as demonstrating carbohydrate-driven community restructuring.

### Block 5 — Results 3.6 / Discussion: dispersion vs location

**SECTION:** Results 3.6 and Discussion (interpretive caveat)
**PRIOR CLAIM:** PERMANOVA rejection can arise from heterogeneous multivariate dispersion rather than a shift in centroid, so homogeneity of dispersion was tested separately and the two results interpreted jointly.
**CITATIONS (lead):** #19 (Anderson & Walsh 2013) · #18 (Anderson 2006) · #20 (Warton et al. 2012) — **(backup):** #22 (Anderson & Santana-Garcon 2015) · #23 (vegan)
**CONFLICTS / NUANCE:** #19 and #20 agree the confounding exists but differ on the remedy (#20 favours model-based multivariate methods over dissimilarity-based tests). Neither test gates the other: a non-significant betadisper does not license an unqualified location interpretation.
**MANUSCRIPT-SAFE SENTENCE:** "Because permutational tests on dissimilarity matrices can respond to differences in multivariate dispersion as well as in centroid location, homogeneity of multivariate dispersions was assessed with betadisper and the two analyses are interpreted together rather than as a sequential decision rule."
**CLAIM-LOCK CHECK:** Supported — and this block is protective, since it constrains rather than expands the compositional claims.

### Block 6 — Methods 2.8: averaging nested technical replicates

**SECTION:** Methods 2.8 (Microbiome data preparation)
**PRIOR CLAIM:** Averaging nested technical/sub-sample replicates to the biological unit before person-level inference is an established and statistically defensible alternative to modelling every sub-sample.
**CITATIONS (lead):** #27 (Lazic 2010) · #26 (Hurlbert 1984) — **(backup):** #28 (Knight et al. 2018) · #30 (Khomich et al. 2021) · #31 (Reichardt et al. 2018, design precedent only)
**CONFLICTS / NUANCE:** Averaging discards information on technical variance, so the magnitude of technical replicate agreement should be reported rather than assumed. **#31 must not be cited for its functional-redundancy conclusion.**
**MANUSCRIPT-SAFE SENTENCE:** "Technical replicate wells were averaged within donor × condition prior to community-level analysis so that inference was conducted at the level of the donor, the independent biological unit of the design."
**CLAIM-LOCK CHECK:** Supported.

### Block 7 — Methods 2.8 / Discussion strengths: resampling sensitivity

**SECTION:** Methods 2.8 and Discussion (strengths)
**PRIOR CLAIM:** Linear mixed models are robust to moderate non-normality and heteroscedasticity of residuals; where diagnostics were imperfect, permutation and bootstrap sensitivity analyses were run to confirm that conclusions did not depend on distributional assumptions.
**CITATIONS (lead):** #32 (Schielzeth et al. 2020) · #34 (Manly 2018) · #33 (Efron & Tibshirani 1993) — **(backup):** #13 (Nakagawa & Schielzeth 2010) · #10 (pbkrtest)
**CONFLICTS / NUANCE:** #32's robustness result applies most strongly to *fixed-effect* estimates and their standard errors; random-effect variance components and their intervals are more sensitive. Because the manuscript reports variance decomposition and ICCs, the robustness claim must not be extended to those quantities — which is the specific reason bootstrap CIs are reported for them.
**MANUSCRIPT-SAFE SENTENCE:** "Because linear mixed models are known to be robust to moderate violations of residual normality and homoscedasticity, the primary models were retained; permutation-based and bootstrap sensitivity analyses were conducted in parallel and yielded concordant conclusions, and bootstrap confidence intervals are reported for all variance components."
**CLAIM-LOCK CHECK:** Supported. **Unsupported if used beyond evidence:** presenting resampling concordance as independent confirmation of a biological effect rather than as robustness of the statistical procedure.

---

## Suggested Methods 2.8 citation block (compact, by theme)

> **LMM.** lme4 [1]; lmerTest [2]; emmeans [3]; nested-data mixed modelling [4, 5]; best practice and pitfalls [6]; variance components [8]; denominator df [9, (10)].
>
> **ICC.** ICC selection and reporting [12]; LMM-based repeatability with bootstrap CIs [13]; rptR [14].
>
> **PERMANOVA.** Method [16]; distance-matrix partitioning and term-specific R² [17]; restricted permutation in nested designs [21]; vegan [23]; Bray–Curtis [24]; dissimilarity theory [25]; effect size/precision [22, 29]; multivariate dispersion [18]; location-vs-dispersion confounding [19, 20].
>
> **Design / replicates.** Experimental unit and pseudoreplication [26, 27]; microbiome design best practice [28]; multivariate method comparison [30].
>
> **Diagnostics.** LMM robustness to distributional violations [32]; permutation [34]; bootstrap [33].

---

## Verification log

- **PubMed-verified (PMID + DOI retrieved):** #6, #12, #13, #18, #28, #31, #15.
- **Crossref-verified (DOI, title, journal, year, authors confirmed):** #1, #2, #3, #4, #5, #7, #8, #9, #10, #11, #14, #16, #17, #19, #20, #21, #22, #23, #24, #25, #26, #27, #29, #30, #32, #33, #34.
- **Not verified — require manual check before use:** Shrout & Fleiss 1979 (*Psychol Bull*); Stoffel, Nakagawa & Schielzeth 2017 (*Methods Ecol Evol* — the rptR journal article as opposed to the CRAN record #14).
- **Deliberately excluded as non-citable:** vegan GitHub issue/discussion threads on `adonis2(by=)` defaults (developer commentary, not literature).
