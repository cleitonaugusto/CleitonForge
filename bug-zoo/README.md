# Bug Zoo — minimal counterexamples for quantum simulator divergences

Public, versioned corpus of minimized circuits on which two
implementations disagree, either two simulator backends or a compiler
pass against the circuit it was given.

Entries come from two places, and the schema differs a little between
them because the questions do:

- The `statevector-*` entries come from [`cforge-fuzz`](../cforge-fuzz),
  the Rust fuzzer, comparing simulator backends. They carry the
  generator seed and the divergence distance at each oracle level.
- The `qiskit-*` entries come from the Python harnesses in
  [`tools/`](../tools), which share the same architecture and run
  against Qiskit's transpiler. They carry the pass name and the version
  it was found on, and the oracle is MQT QCEC rather than an amplitude
  comparison, so there is no distance to record.

Every entry carries the minimal gate list and an OpenQASM 2 reproducer,
which is what makes the finding checkable without rerunning the search.

## Reading an entry

| Field | Meaning |
|---|---|
| `reference_backend` / `device_backend` | The two implementations compared |
| `gates`, `qasm2` | The minimal witness circuit (1-minimal: removing any gate hides the bug) |
| `amplitude_distance` | N1 — max amplitude difference modulo global phase |
| `probability_distance` | N2 — max \|p_a − p_b\| over basis states |
| `classification` | Automatic triage by QGCS dimension |
| `sampling_visible` | `false` = provably invisible to QV, XEB, mirror circuits, RB and any other sampling benchmark, by the conjugation-invariance theorem (see `paper/blindness-hierarchy.tex`) |

Those six fields are the `statevector-*` schema. The `qiskit-*` entries
use `pass`, `qiskit_version`, `oracle` and `criterion` instead, since a
compiler pass is compared against its own input rather than against a
second backend.

One caveat worth stating rather than leaving for someone to discover:
`generator_seed` regenerates the circuit only for entries produced by the
generator version current at the time. `tools/generator.py` had its RNG
stream shifted once, in the commit that added boundary angles, which
broke seed reproduction for `qiskit-16594-single-qubit-003` until it was
restored. `tools/test_generator_stream.py` now guards that. The gate list
and the OpenQASM reproducer never depended on it.

## Fuzzing harnesses

All share one architecture: a weighted random-circuit generator, the N1
amplitude oracle (equality modulo global phase), a greedy 1-minimal
shrinker, and JSON zoo output.

| Harness | Target | Kind |
|---|---|---|
| `cforge-fuzz` (Rust) | native vs quantrs2 / roqoqo / conjugated models | differential |
| `tools/fuzz_qiskit.py` | Qiskit transpiler O3, Aer statevector | metamorphic + differential |
| `tools/fuzz_pytket.py` | tket FullPeepholeOptimise | metamorphic |
| `tools/fuzz_tket_passes.py` | tket RemoveRedundancies / CliffordSimp / SynthesiseTket | metamorphic |
| `tools/fuzz_pennylane.py` | lightning.qubit vs default.qubit | differential |
| `tools/fuzz_pennylane_transforms.py` | cancel_inverses / merge_rotations / single_qubit_fusion / … | metamorphic |

The metamorphic harnesses quotient the transformations the target is
*allowed* to make — global phase always, plus layout permutation and
approximate synthesis for the Qiskit transpiler. A leftover divergence
is then a genuine soundness violation, and the minimal witness names
which quotient (if any) was still missing.

## Regenerating

```bash
cargo run --release -p cforge-fuzz -- \
  --device quantrs2 --oracle amplitude --seed 42 --zoo-dir bug-zoo
python3 tools/fuzz_qiskit.py --iterations 8000 --seed 777 --experiment transpiler
```

## Campaign scoreboard (seed 42, 3000 random circuits each)

| Target | Test | Result |
|---|---|---|
| quantrs2 v0.2.0 | statevector vs native (N1) | **BUG** — Rz sign, 2-gate witness, invisible to all sampling benchmarks |
| Qiskit 2.5.0 transpiler | O3 metamorphic (N1, layout- and approximation-quotiented) | **BUG** — CommutativeCancellation drops sxdg·sxdg = X, 3-gate witness; plus a minor cp small-angle synthesis finding |
| Qiskit 2.5.0 Aer | statevector vs exact (N1) | clean (worst 5.6e-16) |
| tket 2.18.1 | FullPeepholeOptimise metamorphic (N1) | clean (worst 1.0e-13) |
| tket 2.18.1 | RemoveRedundancies / CliffordSimp / SynthesiseTket | clean (16k circuits) |
| PennyLane 0.42 lightning | lightning.qubit vs default.qubit (N1) | clean (worst 3.5e-16) |
| PennyLane 0.42 transforms | cancel_inverses / merge_rotations / single_qubit_fusion / commute_controlled / undo_swaps | clean (20k circuits) |

The one real soundness bug (Qiskit `CommutativeCancellation`) is a
**regression introduced in Qiskit 2.5.0**: releases 1.4.0 through 2.4.0
have `_x_rotations = {"x", "rx"}` (safe); 2.5.0 added `sx` and `sxdg`,
which the cancellation logic mishandles. PennyLane's `cancel_inverses`
and tket's simplification passes — the same class of optimization —
handle the analogous `SX·SX = X` case correctly.

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
