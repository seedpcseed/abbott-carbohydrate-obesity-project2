# Modern-methods environment blockers

**Date:** 2026-07-19

MaAsLin 3 (Bioconductor 3.23, package 1.4.0) and ALDEx3 (CRAN 1.2.0) were prescribed as Priority B analyses in `docs/additional-analyses-plan.2026-07-19.md`.

The render environment used for current Project 2 reports is:

- R 4.3.3
- Bioconductor 3.18
- Maaslin2 1.18.0 / ANCOMBC 2.4.0 installed
- maaslin3 / ALDEx3 **not installed**

Installing BioC 3.23 into this environment would require R ≥4.6 and would break the current knit lock. Per the plan, create a **separate** `renv` project under `integrated/additional_analyses/renv_bioc323/` before running genus-wide qPCR-scaled DA.

Until that environment exists, Absolute-DA outputs remain empty under `absolute_da/`.
