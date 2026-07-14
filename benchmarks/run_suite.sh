#!/usr/bin/env bash
# Canonical benchmark suite — Grover / QFT / BV / QAOA
# Runs each circuit on statevector + quantrs2 and prints comparison table.
set -euo pipefail

CFORGE="${1:-./target/release/cforge}"
CIRCUITS="$(dirname "$0")/circuits"
BACKENDS="statevector,quantrs2"

if ! command -v "$CFORGE" &>/dev/null && [[ ! -x "$CFORGE" ]]; then
    echo "error: cforge binary not found at '$CFORGE'" >&2
    echo "Build with: cargo build --bin cforge --release" >&2
    exit 1
fi

declare -A LABELS=(
    [grover_2q]="Grover (2q)"
    [qft_4q]="QFT (4q)"
    [bv_4q]="BV (4q)"
    [qaoa_maxcut_4q]="QAOA MaxCut (4q)"
)
ORDER=(grover_2q qft_4q bv_4q qaoa_maxcut_4q)

echo ""
echo "=== CleitonForge Canonical Benchmark Suite ==="
echo "  Backend comparison: statevector vs quantrs2"
echo ""

for key in "${ORDER[@]}"; do
    file="$CIRCUITS/${key}.qasm"
    label="${LABELS[$key]}"
    echo "── $label ──────────────────────────"
    "$CFORGE" run --circuit "$file" --backends "$BACKENDS"
    echo ""
done
