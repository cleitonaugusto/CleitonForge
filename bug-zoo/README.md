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

## Campaign scoreboard (seed 42, 3000 random circuits each)

| Target | Test | Result |
|---|---|---|
| quantrs2 v0.2.0 | statevector vs native (N1) | **BUG** — Rz sign, 2-gate witness, invisible to all sampling benchmarks |
| Qiskit 2.5.0 transpiler | O3 metamorphic (N1, layout- and approximation-quotiented) | **BUG** — CommutativeCancellation drops sxdg·sxdg = X, 3-gate witness; plus a minor cp small-angle synthesis finding |
| Qiskit 2.5.0 Aer | statevector vs exact (N1) | clean (worst 5.6e-16) |
| tket 2.18.1 | FullPeepholeOptimise metamorphic (N1) | clean (worst 1.0e-13) |
| PennyLane 0.42 lightning | lightning.qubit vs default.qubit (N1) | see campaign |

The clean results matter as much as the bugs: they show the amplitude
oracle distinguishes a sound optimizer from a broken one, rather than
flagging numerical noise everywhere.

## Current entries

- `statevector-quantrs2-N1-amplitude-001` — the quantrs2 v0.2.0 Rz sign
  inversion, reduced to a 2-gate witness with **zero** probability-level
  signal: the bug class that no standard benchmark can see.
- `statevector-quantrs2-N2-probability-001` — 4-gate witness making the
  same bug visible in measurement counts (bugged rotation mixed with
  correctly-implemented complex gates), rediscovered automatically.
- `statevector-rz-sign-bug-N2-probability-001` — model-backend entry
  demonstrating the strict oracle hierarchy (theorem validation).
- `qiskit-commutative-cancellation-001` — Qiskit 2.5.0 transpiler soundness
  bug: `sxdg·sxdg` (= X) silently cancelled at optimization_level ≥ 2.
  Full standalone reproducer inside the entry.
- `qiskit-cp-small-angle-002` — minor: `cp(θ)` for |θ| ≲ 1e-3 loses both CX
  at O3 even with `approximation_degree=1.0`. A contract/docs question, not a
  wrong-answer bug; kept for completeness.
