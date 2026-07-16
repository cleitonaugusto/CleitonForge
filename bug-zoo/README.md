# Bug Zoo — minimal counterexamples for quantum simulator divergences

Public, versioned corpus of minimized circuits on which two quantum
simulator implementations disagree. Every entry is produced by
[`cforge-fuzz`](../cforge-fuzz) and is fully reproducible: pinned
backend names, the exact gate list, an OpenQASM 2 reproducer, the
generator seed, and the divergence distances at each oracle level.

## Reading an entry

| Field | Meaning |
|---|---|
| `reference_backend` / `device_backend` | The two implementations compared |
| `gates`, `qasm2` | The minimal witness circuit (1-minimal: removing any gate hides the bug) |
| `amplitude_distance` | N1 — max amplitude difference modulo global phase |
| `probability_distance` | N2 — max \|p_a − p_b\| over basis states |
| `classification` | Automatic triage by QGCS dimension |
| `sampling_visible` | `false` = provably invisible to QV, XEB, mirror circuits, RB and any other sampling benchmark, by the conjugation-invariance theorem (see `paper/blindness-hierarchy.tex`) |

## Regenerating

```bash
cargo run --release -p cforge-fuzz -- \
  --device quantrs2 --oracle amplitude --seed 42 --zoo-dir bug-zoo
```

## Current entries

- `statevector-quantrs2-N1-amplitude-001` — the quantrs2 v0.2.0 Rz sign
  inversion, reduced to a 2-gate witness with **zero** probability-level
  signal: the bug class that no standard benchmark can see.
- `statevector-quantrs2-N2-probability-001` — 4-gate witness making the
  same bug visible in measurement counts (bugged rotation mixed with
  correctly-implemented complex gates), rediscovered automatically.
- `statevector-conjugated-*` — model-backend entries demonstrating the
  strict oracle hierarchy (theorem validation).
