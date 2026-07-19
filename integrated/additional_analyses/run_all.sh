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

echo "== Done. Results under $OUT =="
