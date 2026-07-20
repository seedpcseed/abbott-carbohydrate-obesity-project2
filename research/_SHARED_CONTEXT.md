# Shared context for deep-research citation runs

**Project:** Abbott–Lurie ZC07, Project 2 · **Date:** July 19, 2026
**Task:** Assemble leading citations and prior claims to support manuscript drafting (Introduction, Methods, Discussion) **without expanding unsupported claims.** This is authorized academic literature-review support.

## Shared study facts

Human adolescent ex vivo fecal fermentation study (n=16 donors; 8 healthy-weight, 8 with obesity). Each donor's fecal microbiota was incubated with no added carbohydrate, a rapid-digesting carbohydrate formulation (RDC), or a slow-digesting carbohydrate formulation (SDC). Primary endpoints were acetate, propionate, and butyrate net accumulation (0→48 h). 16S V3–V4 community profiling and total bacterial 16S qPCR scaling were secondary/exploratory. This is functional phenotyping of microbial communities, not an in vivo efficacy, digestion, obesity-treatment, or host-glycemic study. Originating rationale: Plaza-Díaz et al. 2022 rodent dietary carbohydrate-replacement study (Front Nutr 2022, DOI 10.3389/fnut.2022.992682).

## Global language guardrails — do NOT retrieve or invent citations that would support:

- SDC superiority over RDC for SCFAs
- obesity equivalence / "preserved fermentation capacity"
- absolute *Fusicatenibacter* expansion
- *Fusicatenibacter* as a proven SDC utilizer or butyrate producer
- median "responder phenotypes"
- microbiome prediction of dietary response
- clinical / glycemic / obesity-treatment benefit of RDC/SDC from this design

## Non-goals for the whole package

These reports are NOT meant to justify: SDC superiority for acetate/propionate/butyrate; obesity-group SCFA equivalence or "preserved capacity"; carbohydrate-driven 48-hour Bray–Curtis restructuring; absolute *Fusicatenibacter* expansion; *Fusicatenibacter* as the mechanistic butyrate source; validated responder phenotypes or dietary-response prediction; clinical/glycemic/obesity-treatment benefit of RDC/SDC.

## Citation integrity rules (MANDATORY)

- Every citation MUST have a real PMID or DOI that you retrieved from a tool (PubMed, Consensus, bioRxiv, or a verifiable web source). **Do not fabricate citations or invent PMIDs/DOIs.**
- If you cannot verify an item, label it explicitly `⚠ UNVERIFIED — needs manual check` rather than presenting it as confirmed.
- Prefer recent reviews (2018–2026) plus seminal primary papers unless the prompt says otherwise.

## Tool usage

Load research tools with ToolSearch, e.g. `ToolSearch("select:mcp__plugin_bio-research_pubmed__search_articles,mcp__plugin_bio-research_pubmed__get_article_metadata")`, and also `mcp__plugin_bio-research_consensus__search`, `mcp__plugin_bio-research_biorxiv__search_preprints`, and `WebSearch` / `WebFetch` as needed. Make multiple targeted searches per report.

## Required output format for every report

```
# Prompt NN — <title>
**Manuscript target:** … · **Citation job:** … · **Date:** 2026-07-19

## Prioritized citation table
| # | Title | Year | Journal | DOI/PMID | Why cite |

## Prior-claim bullets (each mapped to citation numbers)

## Conflicts / nuance

## What we must not overclaim (given claim lock)

## Downstream synthesis blocks
SECTION:
PRIOR CLAIM:
CITATIONS (max 3 lead + 2 backup):
CONFLICTS / NUANCE:
MANUSCRIPT-SAFE SENTENCE:
CLAIM-LOCK CHECK (supported / exploratory / unsupported if used beyond evidence):
```
