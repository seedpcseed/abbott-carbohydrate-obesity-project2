# Prompt 12 — Background fermentation in no-added-carbohydrate controls, and analyte-specific / negative SCFA contrasts

**Manuscript target:** Discussion 4.3 · **Citation job:** Support (a) nontrivial SCFA accumulation in no-added-carbohydrate control wells, (b) mechanistic hypotheses for analyte-specific and negative adjusted contrasts (notably RDC propionate change < control), and (c) claim language that does not force a "more substrate → more SCFA" narrative · **Date:** 2026-07-19

**Citation integrity note:** All 20 entries below were retrieved and verified via the PubMed MCP tool (`search_articles` + `get_article_metadata`) in this session. Every PMID and DOI is copied from the returned metadata record. No item is unverified. Per PubMed tool terms, this report is sourced from PubMed and DOIs are given for every article.

---

## Prioritized citation table

### Tier A — Background / control-well fermentation is expected by assay design

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 1 | An in vitro batch fermentation protocol for studying the contribution of food to gut microbiota composition and functionality | 2021 | Nature Protocols | [10.1038/s41596-021-00537-x](https://doi.org/10.1038/s41596-021-00537-x) · PMID 34089022 | **Strongest single methods citation.** The protocol states explicitly that the oligotrophic medium *with peptone* plus a high concentration of fecal inoculum serve "both to provide the microbiota and **as the main source of nutrients for the bacteria**." Background fermentation is a designed-in feature of batch assays, not an artifact. |
| 2 | The Effect of Amino Acids on Production of SCFA and bCFA by Members of the Porcine Colonic Microbiota | 2022 | Microorganisms | [10.3390/microorganisms10040762](https://doi.org/10.3390/microorganisms10040762) · PMID 35456812 | **Direct quantitative anchor.** A 48-h in vitro batch incubation matched to our timeframe, in which *non-carbohydrate* substrates alone produced ~20 mM acetate, ~10 mM butyrate and ~2 mM propionate. Demonstrates that all three primary analytes accumulate to tens-of-mM without any added carbohydrate. |
| 3 | Dissimilatory amino acid metabolism in human colonic bacteria | 1997 | Anaerobe | [10.1006/anae.1997.0121](https://doi.org/10.1006/anae.1997.0121) · PMID 16887608 | Primary human-fecal batch-culture paper showing acetate, propionate and butyrate are all end products of amino-acid fermentation (e.g. aspartate → acetate + propionate; glutamate → acetate + butyrate; threonine → propionate). Explains peptone-derived background SCFA analyte-by-analyte. |
| 4 | Mucin-Derived O-Glycans Act as Endogenous Fiber and Sustain Mucosal Immune Homeostasis via Short-Chain Fatty Acid Production in Rat Cecum | 2020 | J Nutrition | [10.1093/jn/nxaa097](https://doi.org/10.1093/jn/nxaa097) · PMID 32286621 | Establishes host mucin O-glycans as an **endogenous fermentable substrate** yielding acetate (+37%) and butyrate (+73%). Supports the hypothesis that inoculum-borne host glycans contribute to control-well SCFA. *Rodent in vivo — indirect for our system; label accordingly.* |
| 5 | Enumeration of human colonic bacteria producing phenolic and indolic compounds: effects of pH, carbohydrate availability and retention time on dissimilatory aromatic amino acid metabolism | 1996 | J Appl Bacteriol | [10.1111/j.1365-2672.1996.tb04331.x](https://doi.org/10.1111/j.1365-2672.1996.tb04331.x) · PMID 8810056 | Shows protein/amino-acid fermentation is suppressed ~60% by the presence of a fermentable carbohydrate and ~33% at lower pH. Supplies the *substrate-switching* mechanism by which adding carbohydrate can **redirect** rather than simply add to background flux. |
| 6 | Effect of human faecal donor on in vitro fermentation variables | 1989 | Scand J Gastroenterol | [10.3109/00365528909093060](https://doi.org/10.3109/00365528909093060) · PMID 2544024 | Classic methods citation: donor, substrate, and donor × substrate interactions all significantly affect in vitro SCFA production; recommends ≥3 donors. Justifies our donor-aware modelling and frames donor variance as expected, not anomalous. |

### Tier B — Mechanistic / interpretive: cross-feeding and SCFA utilization

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 7 | Contribution of acetate to butyrate formation by human faecal bacteria | 2004 | Br J Nutrition | [10.1079/BJN20041150](https://doi.org/10.1079/BJN20041150) · PMID 15182395 | **The canonical acetate→butyrate citation.** ¹³C-labelling shows *Faecalibacterium prausnitzii* and *Roseburia* spp. are **net utilisers of acetate**, deriving 85–90% of butyrate carbon from external acetate; 72–91% in mixed faecal *batch* cultures. Directly grounds "SCFAs are consumed as well as produced." |
| 8 | Butyrate production from oligofructose fermentation by the human faecal flora: what is the contribution of extracellular acetate and lactate? | 2006 | Br J Nutrition | PMID 16925864 (no DOI in PubMed record; PII S0007114506002650) | Faecal batch cultures: 80% of newly synthesised butyrate came from **interconversion of extracellular acetate and lactate**, not de novo from the added carbohydrate. The single best citation for "added carbohydrate can shift SCFA interconversion rather than raise all analytes." |
| 9 | Formation of propionate and butyrate by the human colonic microbiota | 2016 | Environ Microbiol | [10.1111/1462-2920.13589](https://doi.org/10.1111/1462-2920.13589) · PMID 27928878 | Definitive mechanistic review of propionate and butyrate pathways, emphasising cross-feeding of lactate, succinate and 1,2-propanediol, and the ecophysiology/environmental responses of the producing taxa. Lead review citation for the whole hypothesis list. |
| 10 | Phylogenetic distribution of three pathways for propionate production within the human gut microbiota | 2014 | ISME J | [10.1038/ismej.2014.14](https://doi.org/10.1038/ismej.2014.14) · PMID 24553467 | Succinate, acrylate and propanediol pathways are held by **different, largely non-overlapping taxa** from butyrate producers — and two Lachnospiraceae (*Coprococcus catus*, *Roseburia inulinivorans*) **switch from butyrate to propionate production on different substrates.** The strongest mechanistic basis for expecting analyte-specific, even opposite-signed, responses. |
| 11 | Lactate-utilizing bacteria, isolated from human feces, that produce butyrate as a major fermentation product | 2004 | Appl Environ Microbiol | [10.1128/AEM.70.10.5810-5817.2004](https://doi.org/10.1128/AEM.70.10.5810-5817.2004) · PMID 15466518 | Key result for substrate competition: **"Addition of glucose to batch cultures prevented lactate utilization until the glucose became exhausted."** Demonstrates that adding a readily fermented carbohydrate can *suspend* a secondary SCFA-generating cross-feeding step — a direct analogue for a lower propionate change under RDC. |
| 12 | Two routes of metabolic cross-feeding between Bifidobacterium adolescentis and butyrate-producing anaerobes from the human gut | 2006 | Appl Environ Microbiol | [10.1128/AEM.72.5.3593-3599.2006](https://doi.org/10.1128/AEM.72.5.3593-3599.2006) · PMID 16672507 | Two distinct cross-feeding mechanisms (consumption of lactate/acetate end products; release of partial breakdown products). Shows carbohydrate structure determines *which* trophic route dominates, hence which SCFA accumulates. |
| 13 | Lactate- and acetate-based cross-feeding interactions between selected strains of lactobacilli, bifidobacteria and colon bacteria in the presence of inulin-type fructans | 2016 | Int J Food Microbiol | [10.1016/j.ijfoodmicro.2016.10.019](https://doi.org/10.1016/j.ijfoodmicro.2016.10.019) · PMID 27810444 | Butyrate formation from a fructan required *acetate* to be present; complete lactate→butyrate conversion occurred only in tricultures. Illustrates that net SCFA output depends on community composition, not substrate supply alone. |
| 14 | Mucin Cross-Feeding of Infant Bifidobacteria and Eubacterium hallii | 2017 | Microb Ecol | [10.1007/s00248-017-1037-4](https://doi.org/10.1007/s00248-017-1037-4) · PMID 28721502 | Explicit statement that **"the ratios of SCFA formed differed depending on the microbial species involved"** in cross-feeding — supports reporting analyte ratios/contrasts separately rather than as a single "fermentation" axis. |
| 15 | The influence of diet on the gut microbiota | 2013 | Pharmacol Res | [10.1016/j.phrs.2012.10.020](https://doi.org/10.1016/j.phrs.2012.10.020) · PMID 23147033 | Accessible review linking SCFA production, cross-feeding on fermentation products, colonic pH lowering, and protein-derived fermentation products. Good general backup citation. |

### Tier C — pH, competition, and community adaptation over the incubation

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 16 | pH and peptide supply can radically alter bacterial populations and short-chain fatty acid ratios within microbial communities from the human colon | 2005 | Appl Environ Microbiol | [10.1128/AEM.71.7.3692-3700.2005](https://doi.org/10.1128/AEM.71.7.3692-3700.2005) · PMID 16000778 | **The pH citation.** Butyrate markedly higher at pH 5.5; **highest propionate at pH 6.5**; and at pH 5.5 high butyrate coincided with *decreasing acetate*. Establishes that a fall in pH driven by fermentation of an added carbohydrate is itself sufficient to **suppress propionate** while favouring butyrate. |
| 17 | The role of pH in determining the species composition of the human colonic microbiota | 2009 | Environ Microbiol | [10.1111/j.1462-2920.2009.01931.x](https://doi.org/10.1111/j.1462-2920.2009.01931.x) · PMID 19397676 | Supplies the taxonomic mechanism for #16: all eight *Bacteroides* spp. tested grew poorly at pH 5.5 (and inhibition worsened in the presence of an SCFA mix), while *E. rectale*-related butyrate producers thrived. Since Bacteroidetes carry the dominant succinate propionate route (#10), acidification predicts a propionate decrement. |
| 18 | Pivotal Roles for pH, Lactate, and Lactate-Utilizing Bacteria in the Stability of a Human Colonic Microbial Ecosystem | 2020 | mSystems | [10.1128/mSystems.00645-20](https://doi.org/10.1128/mSystems.00645-20) · PMID 32900872 | Demonstrates a **genuinely negative** outcome: at pH 5.5, lactate accumulation was concomitant with "greatly reduced butyrate **and propionate** production." Also shows *Bacteroides* and Firmicutes isolates are growth-inhibited by relevant concentrations of lactate **and acetate**. Product-inhibition evidence. |
| 19 | Ecological Adaptation and Succession of Human Fecal Microbial Communities in an Automated In Vitro Fermentation System | 2021 | mSystems | [10.1128/mSystems.00232-21](https://doi.org/10.1128/mSystems.00232-21) · PMID 34313459 | Documents that in vitro communities undergo real compositional turnover: many dominant fecal ASVs fail to maintain dominance while previously undetected "bloomer" ASVs reach significant abundance. Grounds the community-adaptation hypothesis for 24–48 h drift. |
| 20 | Effect of human faecal donor on in vitro fermentation variables | *(see #6)* | | | Cross-listed: donor × substrate interaction. |

### Tier D — Prior examples of non-uniform / absent SCFA responses to added carbohydrate

| # | Title | Year | Journal | DOI / PMID | Why cite |
|---|---|---|---|---|---|
| 21 | Specific substrate-driven changes in human faecal microbiota composition contrast with functional redundancy in short-chain fatty acid production | 2017 | ISME J | [10.1038/ismej.2017.196](https://doi.org/10.1038/ismej.2017.196) · PMID 29192904 | **Closest design analogue.** In vitro *batch* incubations of 15 non-digestible carbohydrates, two pH values, three donors. Some carbohydrates were strongly propionigenic (rhamnose, galactomannans) while fructans and glucans were butyrogenic — i.e. added carbohydrate reallocated SCFAs rather than raising all three. Also reports substantial **functional redundancy**, meaning SCFA output was less substrate-discriminating than composition. |
| 22 | Investigating the response of the butyrate production potential to major fibers in dietary intervention studies | 2024 | npj Biofilms Microbiomes | [10.1038/s41522-024-00533-5](https://doi.org/10.1038/s41522-024-00533-5) · PMID 39080292 | Pooled re-analysis of 14 fiber intervention studies: **"AXOS did not promote an increase in butyrate producers, nor were pathways linked to propionate production stimulated by any intervention."** The single best citation that added carbohydrate routinely fails to move every SCFA — and that propionate is the least responsive. |
| 23 | Understanding the prebiotic potential of different dietary fibers using an in vitro continuous adult fermentation model (PolyFermS) | 2018 | Sci Rep | [10.1038/s41598-018-22438-y](https://doi.org/10.1038/s41598-018-22438-y) · PMID 29531228 | **Irrespective of which fiber was supplemented**, one donor microbiota produced more butyrate and the other more propionate. Donor community, not substrate identity, set the direction of the analyte-specific response. Strong support for donor-aware modelling and for not over-reading a single analyte's contrast. |

---

## Prior-claim bullets (mapped to citation numbers)

- **Batch fermentation media and fecal inocula are themselves fermentable.** Widely used protocols deliberately supply peptone and a high inoculum concentration as the bacteria's main nutrient source [1]; peptone-derived amino acids ferment to acetate, propionate *and* butyrate [3, 15]; and non-carbohydrate substrates alone yield tens of mM of SCFA over a 48-h batch incubation [2].
- **Host-derived glycans in the inoculum are an additional endogenous substrate.** Mucin O-glycans act as "endogenous fiber," raising acetate and butyrate [4]. *(Demonstrated in rodent cecum in vivo, not in our assay.)*
- **Adding a fermentable carbohydrate suppresses background protein fermentation.** Carbohydrate availability reduced dissimilatory amino-acid metabolism by ~60% [5]. Total flux is therefore partly *substituted*, not purely *added*.
- **SCFAs are intermediates, not only end products.** Major butyrate producers are net acetate consumers, deriving 85–90% of butyrate carbon from external acetate [7]; in faecal batch culture 80% of new butyrate came from interconversion of extracellular acetate and lactate [8].
- **Adding a readily fermented sugar can suspend a secondary SCFA-producing step.** Glucose addition prevented lactate utilisation until glucose was exhausted [11] — a documented within-batch catabolite-preference effect.
- **Propionate and butyrate are made by largely distinct taxa via distinct pathways, and some taxa switch product between substrates.** [9, 10]
- **Fermentation-driven acidification favours butyrate and disfavours propionate.** Butyrate is markedly higher and propionate markedly lower at pH 5.5 vs 6.5 [16, 21]; *Bacteroides* — carriers of the dominant succinate propionate route — grow poorly at pH 5.5 [17, 10]. Accumulated lactate and acetate directly inhibit these taxa [18].
- **Genuinely negative SCFA outcomes are documented.** Reduced butyrate *and* propionate under lactate accumulation at low pH [18]; no propionate-pathway stimulation by any of three major fibers across 14 intervention studies [22].
- **Added carbohydrate reallocates rather than uniformly raises SCFAs.** [21, 22, 12, 13, 14]
- **Donor community, not substrate identity, often sets the direction of analyte-specific responses; and SCFA output shows substantial functional redundancy across substrates.** [23, 21, 6]
- **In vitro communities adapt over the incubation.** Dominant fecal taxa lose dominance while low-abundance "bloomers" expand [19].

---

## Graded explanatory hypotheses for our analyte-specific and negative contrasts

Each is a **hypothesis for our data**; the grade describes how directly the cited work demonstrates the *mechanism*, not whether it demonstrates our result.

| # | Hypothesis | Evidence strength | Citations | Grading rationale |
|---|---|---|---|---|
| H1 | Control-well SCFA accumulation reflects fermentation of medium peptone, inoculum-borne residual substrate and host glycans, not assay failure. | **Strong (directly demonstrated in comparable systems)** | 1, 2, 3, 4 | Protocol [1] states the medium/inoculum are the nutrient source by design; [2] quantifies all three analytes from amino acids in a 48-h batch; [3] gives the analyte-level chemistry. Not demonstrated *in our specific medium*. |
| H2 | Adding carbohydrate partly substitutes for background proteolytic fermentation, so net change vs control understates gross carbohydrate fermentation. | **Moderate (mechanism demonstrated; magnitude in our system unknown)** | 5, 3 | [5] directly shows ~60% suppression of protein fermentation by fermentable carbohydrate; the arithmetic consequence for a difference-from-control endpoint is our inference. |
| H3 | Fermentation-driven acidification under an actively fermented carbohydrate selectively suppresses propionate formation (via inhibition of *Bacteroides*/succinate-pathway taxa) while sparing or favouring butyrate. | **Moderate–strong (well demonstrated in continuous and batch culture; pH not manipulated here)** | 16, 17, 10, 21 | [16] and [21] show the propionate/butyrate pH split experimentally; [17] and [10] supply the taxonomic and pathway mechanism. **We did not measure or control pH**, so this remains a hypothesis for our assay. |
| H4 | Net propionate is reduced because propionate and its precursors (succinate, lactate, 1,2-propanediol) are consumed or rerouted, and/or because product inhibition by accumulating acetate and lactate constrains propionigenic taxa. | **Moderate (cross-feeding and product inhibition demonstrated; net consumption not measured here)** | 9, 10, 18, 11 | [18] documents concurrent butyrate and propionate reduction plus acetate/lactate growth inhibition; [9, 10] establish the precursor cross-feeding routes. We have no isotope or intermediate data. |
| H5 | Rapid availability of a readily fermented carbohydrate transiently suspends secondary cross-feeding steps (catabolite preference), lowering the SCFAs that depend on those steps. | **Moderate (directly demonstrated in batch culture for lactate utilisation)** | 11, 13 | [11] is an explicit batch-culture demonstration of the phenomenon, though for lactate→butyrate rather than propionate. Extension to propionate is inference. |
| H6 | Because acetate is simultaneously the largest product and a major consumed substrate, acetate change is a poor index of total fermentation and can move independently of butyrate. | **Strong (isotopically demonstrated)** | 7, 8, 16 | ¹³C tracer work [7, 8] establishes net acetate uptake and interconversion; [16] shows acetate falling as butyrate rises. |
| H7 | Producer taxa are partly interchangeable, so community-level SCFA output is functionally redundant and less discriminating between carbohydrates than composition is — reducing power to detect active-vs-control differences. | **Moderate–strong (demonstrated across donors and substrates)** | 21, 23, 22 | [21] reports this explicitly; [23] shows donor sets analyte direction irrespective of fiber; [22] shows null propionate responses across 14 studies. |
| H8 | Community composition drifts over 24–48 h, so the 48-h endpoint reflects an adapted community whose product profile need not track initial substrate supply. | **Weak–moderate (demonstrated in longer continuous systems; extrapolated to 48 h batch)** | 19, 21 | [19] used 22–31-day continuous fermentors, so the timescale is a substantial extrapolation. Flag as speculative. |
| H9 | Between-donor heterogeneity is a dominant variance component, so modest fixed carbohydrate effects are hard to resolve. | **Strong (long-established methodological finding)** | 6, 23, 21 | [6] established donor and donor × substrate effects in 1989 and recommended multiple donors; [23] and [21] confirm in modern systems. |

---

## Conflicts / nuance

- **Direction of the pH effect is consistent, but our pH was not manipulated.** [16] and [21] both report more butyrate at low pH and more propionate at high pH. However, both *set* pH as an experimental factor. In an unbuffered or weakly buffered batch assay, pH is a *consequence* of fermentation; we did not record it. H3 must not be stated as an established cause of our propionate result.
- **"Functional redundancy" cuts both ways.** [21] found SCFA output "surprisingly reproducible" across donors — which argues that a real substrate effect *should* have been detectable. That is a genuine tension with a purely power-based reading of our null contrasts, and should be acknowledged rather than suppressed.
- **Species-level attribution is unsafe.** [10] shows *Coprococcus catus* and *Roseburia inulinivorans* switch between butyrate and propionate depending on substrate. This is a reason to avoid mapping any observed SCFA change onto a named taxon in our 16S data.
- **Timescale mismatch.** Nearly all cross-feeding evidence [7, 8, 11, 12, 13] comes from pure/co-culture or short batch systems; the adaptation evidence [19] comes from multi-week fermentors. Neither is a perfect match to a 48-h complex-community batch.
- **Cross-species evidence.** [2] is porcine-derived microbiota and [4] is rodent in vivo. Both are mechanistically informative but must be labelled as non-human when cited.
- **A single non-significant analyte-specific negative contrast is fragile.** The RDC propionate decrement is one contrast among several tested; it should be presented as hypothesis-generating and consistent with known mechanisms, not as a demonstrated suppression effect.

---

## What we must not overclaim (given claim lock)

- ❌ **Do not** state that SDC and RDC produce equivalent SCFA responses, or that they are "comparable." Non-significant SDC-vs-RDC contrasts mean the two were **not distinguished under this design and sample size** — an absence of evidence, not evidence of equivalence. No equivalence test was performed.
- ❌ **Do not** claim SDC superiority over RDC for any SCFA. Nothing here supports it and the prompt's own results do not.
- ❌ **Do not** convert the RDC propionate decrement into a claim that RDC "suppresses" or "reduces" propionate production. It is a lower *net change vs control* in adjusted contrasts, in a design where the control itself ferments actively.
- ❌ **Do not** attribute the analyte-specific pattern to pH. We did not measure pH. H3 is a hypothesis with external support only.
- ❌ **Do not** attribute any SCFA change to *Fusicatenibacter* or to any named taxon as a utilizer or producer.
- ❌ **Do not** frame background control-well fermentation as a limitation that undermines the assay; the literature treats it as an expected property. Equally, do not use it to explain away every null result.
- ❌ **Do not** extend any of this to obesity-group comparisons, host glycemic outcomes, digestion kinetics in vivo, or clinical benefit.
- ⚠️ **Label as exploratory:** all of H1–H9 when applied to our data; the community-adaptation hypothesis (H8) is the weakest and should be presented as speculative.

---

## Discussion 4.3 paragraph scaffold

> **[Sentence 1 — state the observation plainly, without a deficit framing.]**
> Acetate, propionate and butyrate accumulated over 48 h in every condition, including wells to which no carbohydrate was added.
>
> **[Sentence 2 — normalize it against assay design. Cite 1, 2, 3.]**
> This is an expected property of fecal batch-fermentation assays rather than an indication of assay failure: widely used protocols supply a peptone-containing medium and a high-density fecal inoculum that together serve as the principal nutrient source for the community [1], and non-carbohydrate substrates alone are sufficient to generate tens of millimolar concentrations of all three analytes over a comparable 48-h incubation [2], with amino-acid fermentation yielding acetate, propionate and butyrate in substrate-dependent proportions [3].
>
> **[Sentence 3 — state the arithmetic consequence for our endpoint. Cite 5.]**
> Because a fermentable carbohydrate also suppresses background proteolytic fermentation [5], the change-from-control endpoint used here reflects the *net* difference between two actively fermenting systems, and may understate gross carbohydrate-derived flux.
>
> **[Sentence 4 — introduce the analyte-specific finding without editorializing.]**
> Against this background, active-versus-control contrasts were largely non-significant, and the adjusted propionate change under RDC was lower than under no added carbohydrate.
>
> **[Sentence 5 — the core interpretive move: SCFAs are intermediates. Cite 7, 8, 9.]**
> Interpreting such contrasts requires recognizing that short-chain fatty acids are metabolic intermediates as well as end products: major butyrate producers are net consumers of acetate, deriving the large majority of butyrate carbon from extracellular acetate [7], and in fecal batch culture most newly formed butyrate arises from interconversion of extracellular acetate and lactate rather than de novo from the added substrate [8]. Propionate in particular is generated by taxa and pathways largely distinct from those producing butyrate [9, 10].
>
> **[Sentence 6 — offer the hypotheses, explicitly labelled. Cite 16, 17, 11, 18.]**
> Several non-exclusive mechanisms could contribute to an analyte-specific or negative propionate contrast, though none can be adjudicated with the present data: fermentation-driven acidification, which in controlled systems raises butyrate while lowering propionate and inhibits the *Bacteroides* spp. that carry the dominant succinate route [16, 17]; transient suspension of secondary cross-feeding steps in the presence of a readily fermented sugar [11]; and inhibition of propionigenic taxa by accumulating acetate and lactate [18]. We did not measure pH or fermentation intermediates, so these remain hypotheses.
>
> **[Sentence 7 — situate within prior null/non-uniform findings. Cite 22, 21, 23.]**
> Non-uniform SCFA responses to added carbohydrate are well precedented: a pooled analysis of fourteen fiber-intervention studies found no stimulation of propionate-associated pathways by any of three major fibers [22], and in vitro work across multiple carbohydrates and donors shows substrates reallocating SCFA profiles — propionigenic versus butyrogenic — rather than raising all analytes together [21], with the donor community often determining the direction of response irrespective of which fiber is supplied [23].
>
> **[Sentence 8 — the SDC/RDC statement, worded to avoid equivalence. Cite 21, 6.]**
> Differences between SDC and RDC were not statistically resolved for any of the three primary short-chain fatty acids. Given the marked between-donor variation characteristic of this assay [6] and the functional redundancy in short-chain fatty acid output reported across donors and substrates [21], we interpret this as a limit of resolution under the present design rather than as evidence that the two formulations are metabolically equivalent; no equivalence testing was performed.
>
> **[Sentence 9 — close on scope.]**
> These are ex vivo functional measurements of microbial communities and do not speak to in vivo digestion, host glycemic response, or clinical outcome.

---

## Downstream synthesis blocks

**SECTION:** Discussion 4.3 — background fermentation in no-added-carbohydrate controls
**PRIOR CLAIM:** Fecal batch-fermentation systems generate substantial SCFA in the absence of added test carbohydrate, because medium components (notably peptone) and the fecal inoculum itself are fermentable.
**CITATIONS (lead):** Pérez-Burillo 2021 Nat Protoc (PMID 34089022); Van den Abbeele 2022 Microorganisms (PMID 35456812); Smith & Macfarlane 1997 Anaerobe (PMID 16887608). **(backup):** Hino 2020 J Nutr (PMID 32286621); Smith & Macfarlane 1996 J Appl Bacteriol (PMID 8810056).
**CONFLICTS / NUANCE:** The quantitative anchor [PMID 35456812] used porcine-derived microbiota; the mucin evidence [PMID 32286621] is rodent in vivo. Neither replicates our medium. Magnitude of background flux in our specific system is not established by these citations.
**MANUSCRIPT-SAFE SENTENCE:** "SCFA accumulation in no-added-carbohydrate wells is an expected feature of fecal batch fermentation, in which the peptone-containing medium and the fecal inoculum themselves serve as fermentable substrate."
**CLAIM-LOCK CHECK:** Supported.

---

**SECTION:** Discussion 4.3 — why net contrasts understate gross fermentation
**PRIOR CLAIM:** Added fermentable carbohydrate suppresses background proteolytic fermentation, so a difference-from-control endpoint captures net rather than gross substrate-derived flux.
**CITATIONS (lead):** Smith & Macfarlane 1996 J Appl Bacteriol (PMID 8810056); Smith & Macfarlane 1997 Anaerobe (PMID 16887608). **(backup):** Scott 2013 Pharmacol Res (PMID 23147033).
**CONFLICTS / NUANCE:** The ~60% suppression figure was measured for phenolic/indolic products of aromatic amino acids, not for total SCFA. The extension to our endpoint arithmetic is our inference.
**MANUSCRIPT-SAFE SENTENCE:** "Because fermentable carbohydrate also reduces background amino-acid fermentation, change-from-control endpoints reflect the net difference between two actively fermenting systems."
**CLAIM-LOCK CHECK:** Exploratory (mechanism supported; quantitative extension to our endpoint is inference).

---

**SECTION:** Discussion 4.3 — SCFAs as intermediates, not only end products
**PRIOR CLAIM:** Acetate and lactate are consumed and interconverted to butyrate, so individual SCFA changes do not index total fermentation.
**CITATIONS (lead):** Duncan 2004 Br J Nutr (PMID 15182395); Morrison 2006 Br J Nutr (PMID 16925864); Louis & Flint 2016 Environ Microbiol (PMID 27928878). **(backup):** Belenguer 2006 AEM (PMID 16672507); Moens 2016 Int J Food Microbiol (PMID 27810444).
**CONFLICTS / NUANCE:** Isotope evidence comes from pure/co-culture and short batch systems; we have no tracer or intermediate measurements to confirm interconversion in our incubations.
**MANUSCRIPT-SAFE SENTENCE:** "Short-chain fatty acids are metabolic intermediates as well as end products; butyrate-producing taxa are net consumers of acetate, and most newly formed butyrate in fecal batch culture derives from interconversion of extracellular acetate and lactate."
**CLAIM-LOCK CHECK:** Supported (as a general statement about the system); exploratory if applied to explain a specific contrast in our data.

---

**SECTION:** Discussion 4.3 — mechanisms for the lower RDC propionate change
**PRIOR CLAIM:** Analyte-specific or negative propionate contrasts are mechanistically plausible via acidification-driven suppression of succinate-pathway taxa, catabolite-preference suspension of cross-feeding, and product inhibition.
**CITATIONS (lead):** Walker 2005 AEM (PMID 16000778); Duncan 2009 Environ Microbiol (PMID 19397676); Wang 2020 mSystems (PMID 32900872). **(backup):** Duncan 2004 AEM (PMID 15466518); Reichardt 2014 ISME J (PMID 24553467).
**CONFLICTS / NUANCE:** pH was an experimentally *set* variable in the cited work but was neither controlled nor measured in ours. The catabolite-preference demonstration [PMID 15466518] concerns lactate→butyrate, not propionate. All three mechanisms are untested here and are not mutually exclusive.
**MANUSCRIPT-SAFE SENTENCE:** "A lower net propionate change under RDC is consistent with several established mechanisms — fermentation-driven acidification disfavouring succinate-pathway propionate producers, transient suspension of secondary cross-feeding, and product inhibition by accumulating acetate and lactate — but we did not measure pH or fermentation intermediates and cannot distinguish among them."
**CLAIM-LOCK CHECK:** Exploratory. **Unsupported if used** to assert that RDC causally suppresses propionate production, or that pH mediated the effect.

---

**SECTION:** Discussion 4.3 — precedent for non-uniform SCFA responses
**PRIOR CLAIM:** Added carbohydrate frequently fails to raise all SCFAs uniformly; propionate is the least responsive analyte, and donor community often sets response direction.
**CITATIONS (lead):** Van-Wehle & Vital 2024 npj Biofilms Microbiomes (PMID 39080292); Reichardt 2017 ISME J (PMID 29192904); Poeker 2018 Sci Rep (PMID 29531228).
**CONFLICTS / NUANCE:** [PMID 29192904] also reports SCFA production as "surprisingly reproducible" across donors, indicating real substrate effects are often detectable — a tension with reading our nulls purely as low power. [PMID 39080292] pools in vivo intervention studies, not batch fermentations.
**MANUSCRIPT-SAFE SENTENCE:** "Non-uniform SCFA responses are well precedented: propionate-associated pathways were not stimulated by any of three major fibers across fourteen pooled intervention studies, and in vitro substrates commonly reallocate the SCFA profile rather than raising all analytes together."
**CLAIM-LOCK CHECK:** Supported.

---

**SECTION:** Discussion 4.3 — SDC versus RDC
**PRIOR CLAIM:** Non-significant SDC-vs-RDC contrasts reflect the resolution limit of this design, not demonstrated equivalence.
**CITATIONS (lead):** McBurney & Thompson 1989 Scand J Gastroenterol (PMID 2544024); Reichardt 2017 ISME J (PMID 29192904). **(backup):** Poeker 2018 Sci Rep (PMID 29531228).
**CONFLICTS / NUANCE:** Functional redundancy [PMID 29192904] argues both for expecting reproducible SCFA output *and* for expecting small between-substrate differences — it does not by itself license an equivalence reading. No equivalence testing was performed in our study.
**MANUSCRIPT-SAFE SENTENCE:** "Differences between SDC and RDC were not statistically resolved for any primary short-chain fatty acid; given documented between-donor variation and functional redundancy in SCFA output, we interpret this as a limit of resolution under the present design rather than as evidence of metabolic equivalence."
**CLAIM-LOCK CHECK:** Supported as worded. **Unsupported if used** to claim equivalence, "preserved fermentation capacity," or SDC superiority.

---

**SECTION:** Discussion 4.3 — community adaptation over the incubation
**PRIOR CLAIM:** In vitro fecal communities restructure over the incubation, so a 48-h endpoint reflects an adapted community.
**CITATIONS (lead):** Gnanasekaran 2021 mSystems (PMID 34313459). **(backup):** Reichardt 2017 ISME J (PMID 29192904).
**CONFLICTS / NUANCE:** The adaptation evidence comes from 22–31-day continuous fermentors; extrapolation to a 48-h batch incubation is substantial and unvalidated.
**MANUSCRIPT-SAFE SENTENCE:** "In vitro fecal communities are known to restructure during incubation, with some dominant fecal taxa declining and low-abundance taxa expanding; the extent to which this occurs within 48 h in our system is unknown."
**CLAIM-LOCK CHECK:** Exploratory / speculative — weakest hypothesis in this set; include only with explicit hedging, or omit.
