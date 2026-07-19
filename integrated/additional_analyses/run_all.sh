#!/usr/bin/env bash
# Run Priority A/B additional analyses that are possible in the current environments.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export ABBOTT_ROOT="$ROOT"
OUT="$ROOT/integrated/results/additional_analyses_2026-07-19"
mkdir -p "$OUT/diagnostics"

echo "== Phase A: matched table =="
python3 "$ROOT/integrated/additional_analyses/python/build_matched_table.py"

echo "== Phase B: focal Fusicatenibacter + total load =="
python3 "$ROOT/integrated/additional_analyses/python/run_focal_and_load.py"

echo "== Phase: SCFA contrasts, alpha, PERMANOVA =="
RBIN="${ABBOTT_RBIN:-/tmp/abbott-render-r/bin/Rscript}"
if [[ ! -x "$RBIN" ]]; then
  RBIN="$(command -v Rscript)"
fi
"$RBIN" "$ROOT/integrated/additional_analyses/R/run_scfa_alpha_permanova.R" \
  2>&1 | tee "$OUT/diagnostics/r_scfa_alpha_permanova.log"

echo "== Phase: genus matrices + MaAsLin3 / ALDEx3 (R>=4.6 BioC 3.23 lib) =="
python3 "$ROOT/integrated/additional_analyses/python/prepare_genus_matrices.py"
BIOC_LIB="$ROOT/integrated/additional_analyses/renv_bioc323/library"
if [[ -d "$BIOC_LIB" ]] && command -v Rscript >/dev/null; then
  export R_LIBS="$BIOC_LIB"
  export R_LIBS_USER="$BIOC_LIB"
  Rscript "$ROOT/integrated/additional_analyses/R/absolute_da_maaslin3.R" "$ROOT" \
    2>&1 | tee "$OUT/absolute_da/maaslin3_run.log"
  Rscript "$ROOT/integrated/additional_analyses/R/absolute_da_aldex3.R" "$ROOT" \
    2>&1 | tee "$OUT/absolute_da/aldex3_run.log"
else
  echo "Skipping MaAsLin3/ALDEx3: missing $BIOC_LIB or system Rscript" | tee "$OUT/absolute_da/SKIPPED.md"
fi

echo "== Done. Results under $OUT =="
