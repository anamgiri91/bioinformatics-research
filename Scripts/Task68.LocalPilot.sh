#!/bin/bash
# Local, sequential execution. No SSH, remote scheduler or cluster access.
# Usage: bash Scripts/Task68.LocalPilot.sh ROW [--dry-run]
set -euo pipefail
SRN_LOCAL_REPO="$(cd "$(dirname "$0")/.." && pwd)"
SRN_LOCAL_ROW="${1:-1}"
SRN_LOCAL_MODE="${2:-}"
case "$SRN_LOCAL_ROW" in 1|2) ;; *) echo 'Pilot row must be 1 or 2' >&2; exit 2;; esac
case "$SRN_LOCAL_MODE" in ''|--dry-run) ;; *) echo 'Only --dry-run is supported as the second argument' >&2; exit 2;; esac
export PATH="$SRN_LOCAL_REPO/.venv-shared-fragment/bin:$SRN_LOCAL_REPO/Data/shared_fragment_tools/Bismark-0.24.2:$SRN_LOCAL_REPO/Data/shared_fragment_tools/TrimGalore:$PATH"
cd "$SRN_LOCAL_REPO"
SRN_LOCAL_ARGS=(
  --manifest Results/Task68_SpermPilotSamples.tsv --row "$SRN_LOCAL_ROW"
  --work Data/gse165915_sperm_wgbs
  --reference-dir Data/gse165915_sperm_reference/hg38_bismark
  --map-dir Data/gse165915_sperm_reference/hg38_cpg_map
  --threads 4
)
if [ "$SRN_LOCAL_MODE" = '--dry-run' ]; then SRN_LOCAL_ARGS+=(--dry-run); fi
exec "$SRN_LOCAL_REPO/.venv-shared-fragment/bin/python" Scripts/Task68.SpermPilot.py "${SRN_LOCAL_ARGS[@]}"
