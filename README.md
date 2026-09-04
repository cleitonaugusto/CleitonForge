# CleitonForge

**[English]** | [Português](README.pt-BR.md)

[![GitHub Sponsors](https://img.shields.io/github/sponsors/cleitonaugusto?label=Sponsor&logo=GitHub&color=ea4aaa)](https://github.com/sponsors/cleitonaugusto)
[![Crates.io](https://img.shields.io/crates/v/cleitonforge)](https://crates.io/crates/cleitonforge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22307398.svg)](https://doi.org/10.5281/zenodo.22307398)
[![PyPI](https://img.shields.io/pypi/v/cleitonforge)](https://pypi.org/project/cleitonforge/)

**CleitonForge** is a differential fuzzer for quantum compilers and simulators,
written in Rust. It generates random circuits, runs each one through more than
one implementation, and reports the cases where the implementations disagree.
When it finds one, it shrinks the circuit down to the smallest version that
still fails. The CLI binary is `cforge`.

It found a soundness bug in Qiskit's transpiler:
[issue #16594](https://github.com/Qiskit/qiskit/issues/16594). The
`CommutativeCancellation` pass cancelled `sxdg sxdg sx`, which is an X, down to
an empty circuit at optimization level 2 and above. No error, no warning, and
the measured result changes. A core maintainer confirmed it and it was fixed in
Qiskit 2.5.1.

---

## The problem

Simulators and compilers are hard to test because there is usually nothing to
compare an answer against. A wrong amplitude is still a well formed amplitude.
Nothing raises, nothing logs, and the number you get back looks exactly like a
correct one. You only notice when you already know what the answer should be,
which is the case you did not need the tool for.

Differential testing gets around this. Instead of asking "is this right", it
asks "do these two agree". Two implementations of the same specification should
return the same state for the same circuit. When they do not, one of them is
wrong, and the disagreement is the signal.

That is the whole idea. The rest is making it precise enough to be worth
reading:

- **The generator** aims at the numeric thresholds a compiler branches on:
  multiples of pi/2, angles just under a synthesis cutoff, sequences that
  collapse to identity. Uniform random angles almost never land on those.
- **The oracle** compares modulo global phase, because global phase is not
  observable and flagging it would bury the real findings under noise. The
  harnesses for external compilers can hand the verdict to
  [MQT QCEC](https://github.com/cda-tum/mqt-qcec) instead, which reasons
  symbolically rather than simulating.
- **The shrinker** reduces a failure to a 1-minimal witness, meaning removing
  any single gate makes the bug disappear. #16594 came out of it as three gates
  on one qubit.
- **The triage** cross-checks every verdict against the exact operator before
  calling anything a bug, and tells a known class from a new one by accumulated
  rotation instead of by gate name.

## The hit rate, honestly

Roughly one real bug per hundred thousand circuits. These stacks are mature and
well tested, and most disagreements turn out to be a declared transformation
that the oracle had not been told about, not a fault. Campaigns against tket,
PennyLane and Cirq came back clean.

So the value here is not volume. It is that the findings which survive triage
are real, minimal, and reproducible by someone who does not trust you. They live
in [`bug-zoo/`](bug-zoo), one JSON file each, with the gate list, an OpenQASM 2
reproducer and a record of how it was found. The backend divergences also carry
the generator seed and the distance at each oracle level.

## Why Rust

- The generator, the state vector and the shrinker are all in the hot path of a
  campaign that runs tens of thousands of circuits.
- Qiskit's transpiler is itself moving to Rust, and the bug above is in Rust
  code. It helps to read the language of the thing you are fuzzing.
- For the benchmarking side, predictable low overhead matters, because
  measurement noise is what is being measured.

## The benchmarking side

CleitonForge started as a cross-backend benchmarking layer, and that part still
works and is documented below. `cforge run` executes the same OpenQASM file on
more than one backend and reports time, memory, depth and fidelity side by side.
It is now a component of the fuzzer rather than the point of the project.

Two rules kept from that period, because they are what makes the comparison
mean anything:

- Never optimized to favor a specific simulator.
- Plugin architecture (`SimulationBackend` trait), so a backend can be added
  without touching the core IR.

---

## Architecture

```
CleitonForge/                    Rust workspace
├── cforge-core/                 Canonical IR: Circuit, GateKind, Operation
├── cforge-parser/               OpenQASM 2 + OpenQASM 3 parsers
├── cforge-backends/             SimulationBackend trait + implementations
├── cforge-metrics/              Fidelity, depth, timing, memory measurement
├── cforge-fuzz/                 The fuzzer
│   ├── generator.rs             Weighted circuit generator
│   ├── oracle.rs                Divergence oracles (N1 amplitude, N2 probability)
│   ├── shrinker.rs              Greedy reduction to a 1-minimal witness
│   ├── triage.rs                Known class vs new finding
│   └── zoo.rs                   JSON output for bug-zoo/
├── cforge-cli/                  `cforge` binary (clap + comfy-table)
│   └── examples/
│       └── compare_grover.rs   Grover algorithm, pure Rust API example
├── cforge-py/                   Python bindings (PyO3)
├── tools/                       Harnesses for external compilers
│   ├── fuzz_qiskit*.py          Qiskit: full pipeline, isolated passes, wide
│   ├── fuzz_pytket*.py          tket
│   ├── fuzz_pennylane*.py       PennyLane devices and transforms
│   ├── fuzz_cirq.py             Cirq
│   ├── oracle_qcec.py           MQT QCEC as a three way verdict
│   └── triage_known.py          Triage by accumulated rotation
├── bug-zoo/                     Minimized counterexamples, one JSON each
└── examples/
    └── bell.qasm               Bell state in OpenQASM 2
```

### Data flow

```
.qasm file
    │
    ▼
cforge-parser  ──►  Circuit (canonical IR)
                         │
              ┌──────────┴──────────┐
              ▼                     ▼
  NativeStateVector          QuantRS2Backend
     Backend                       │
              └──────────┬──────────┘
                         ▼
              cforge-metrics: fidelity, time, memory
                         │
                         ▼
              cforge-cli: table / JSON output
```

---

## Quick Start

### Prerequisites

- Rust 1.96+ (`rustup update stable`)
- No external simulators required. All dependencies are pure Rust crates.

### Build

```bash
git clone https://github.com/cleitonaugusto/CleitonForge.git
cd CleitonForge
cargo build --release
# binary at: target/release/cforge
```

### Run a circuit on both backends

```bash
cforge run --circuit examples/bell.qasm --backends statevector,quantrs2 --shots 1024
```

```
Circuit: 2 qubits  |  2 gates  |  depth 2  |  seed 0xdeadbeefcafebabe
┌──────────────────────┬───────────┬────────┬───────┬───────┬──────────┬───────┐
│ Backend              │ Time (ms) │ Memory │ Depth │ Gates │ Fidelity │ Shots │
╞══════════════════════╪═══════════╪════════╪═══════╪═══════╪══════════╪═══════╡
│ statevector-native   │ 0.002     │ 64 B   │ 2     │ 2     │ 1.000000 │ 1024  │
├──────────────────────┼───────────┼────────┼───────┼───────┼──────────┼───────┤
│ statevector-quantrs2 │ 0.004     │ 64 B   │ 2     │ 2     │ 1.000000 │ 1024  │
└──────────────────────┴───────────┴────────┴───────┴───────┴──────────┴───────┘
```

### Validate a circuit without simulation

```bash
cforge validate --circuit examples/bell.qasm
```

```
File   : examples/bell.qasm
Qubits : 2
Gates  : 2
Depth  : 2
By gate:
  cx       1
  h        1
Status : OK
```

### Export results as JSON (for scripts / CI)

```bash
cforge run --circuit examples/bell.qasm --shots 1024 --format json
```

```json
{
  "circuit": { "file": "examples/bell.qasm", "qubits": 2, "gates": 2, "depth": 2 },
  "shots": 1024,
  "seed": 16935479246996842942,
  "results": [
    {
      "backend": "statevector-native",
      "time_ms": 0.456,
      "memory_bytes": 64,
      "depth": 2,
      "gates": 2,
      "fidelity": 1.0,
      "shots": 1024,
      "error": null
    }
  ]
}
```

### Grover's algorithm end-to-end example

3-qubit Grover search for |101⟩, 2 iterations, both backends:

```bash
cargo run --example compare_grover -p cforge-cli
```

```
Target state : |101⟩  (index 5)
Circuit      : 43 gates  |  depth 21

Backend  : statevector-native
  Top state  : |101⟩  prob = 0.9453 (94.5 %)   [theory: sin²(5θ) ≈ 94.8 %]
  Fidelity   : 1.00000000

Cross-backend fidelity (native vs quantrs2): 1.00000000
Both backends agree: YES ✓
```

---

## CLI Reference

```
cforge run
  --circuit <path>          OpenQASM 2 or 3 file (auto-detected)
  --backends <list>         Comma-separated: statevector, quantrs2  [default: both]
  --shots <n>               Measurement shots; 0 = statevector only  [default: 0]
  --seed <u64>              PRNG seed for reproducible counts  [default: 0xdeadbeef_cafebabe]
  --format <table|json>     Output format  [default: table]

cforge validate
  --circuit <path>          Parse and report circuit statistics; exit 1 if invalid
```

---

## Supported Gate Set

CleitonForge implements the full **OpenQASM 3 `stdgates.inc`** gate set:

| Category | Gates |
|---|---|
| Single-qubit, no params | `id` `x` `y` `z` `h` `s` `sdg` `t` `tdg` `sx` `sxdg` |
| Single-qubit, parametric | `rx(θ)` `ry(θ)` `rz(θ)` `p(θ)` `u(θ,φ,λ)` |
| Two-qubit | `cx` `cy` `cz` `ch` `csx` `crx` `cry` `crz` `cp` `cu` `swap` |
| Three-qubit | `ccx` (Toffoli) `cswap` (Fredkin) |
| Aliases | `cnot`→`cx` `u1`→`p` `u2(φ,λ)`→`u(π/2,φ,λ)` `u3`→`u` `ccnot`→`ccx` `fredkin`→`cswap` |

---

## Supported Input Formats

| Format | Auto-detected by | Notes |
|---|---|---|
| OpenQASM 2.0 | everything else | `include` resolved from file directory; no stdgates disk file needed |
| OpenQASM 3 | `OPENQASM 3` header | Parses `stdgates.inc` gate calls via `oq3_semantics` |

Both parsers support **whole-register gate application**:
`h q;` on a 3-qubit register expands to `h q[0]; h q[1]; h q[2];`

---

## Rust API

```rust
use cforge_core::{Circuit, GateKind, Operation};
use cforge_backends::{DEFAULT_SEED, NativeStateVectorBackend, SimulationBackend};
use cforge_metrics::{compute_stats, measure};

// Build a Bell state circuit
let mut circuit = Circuit::new(2);
circuit.push(Operation::new(GateKind::H,  vec![0], vec![]));
circuit.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));

// Run and measure
let stats = compute_stats(&circuit);
let result = NativeStateVectorBackend.run(&circuit, 1024, DEFAULT_SEED)?;

println!("depth = {}", stats.depth);               // 2
println!("|00⟩ prob = {:.3}", result.statevector[0].norm_sqr()); // 0.500
```

---

## Adding a New Backend

Implement the `SimulationBackend` trait (one method):

```rust
use cforge_backends::{BackendError, SimulationBackend, SimulationResult};
use cforge_core::Circuit;

pub struct MyBackend;

impl SimulationBackend for MyBackend {
    fn name(&self) -> &str { "my-backend" }

    fn run(
        &self,
        circuit: &Circuit,
        shots: usize,
        seed: u64,
    ) -> Result<SimulationResult, BackendError> {
        // ... your simulation logic
    }
}
```

Then pass `--backends my-backend` to the CLI after registering it in `cforge-cli/src/main.rs`.

---

## Memory Measurement

For statevector simulations the peak memory is:

```
2^n_qubits × 16 bytes  (two f64 per Complex128 amplitude)
```

On Linux, CleitonForge measures the actual RSS delta via `/proc/self/status`
while the statevector is live. For circuits smaller than one OS page (4 KiB)
the theoretical value is used instead. Both values are included in JSON output.

| Qubits | Theoretical peak |
|--------|-----------------|
| 10     | 16 KB           |
| 20     | 16 MB           |
| 22     | 64 MB (max)     |

---

## Status

| Phase | Crate             | Status |
|-------|-------------------|--------|
| 0     | workspace setup   | ✓      |
| 1     | `cforge-core`     | ✓      |
| 2     | `cforge-parser`   | ✓      |
| 3     | `cforge-backends` | ✓      |
| 4     | `cforge-metrics`  | ✓      |
| 5     | `cforge-cli`      | ✓      |
| 6     | examples + docs   | ✓      |
| 7     | `cforge-py`       | ✓      |
| 8     | `cforge-fuzz`     | ✓      |

Python bindings are published: `pip install cleitonforge`.

**Planned:** additional backends (qoqo, q1tsim), noise modeling, extended
OpenQASM 3 gate coverage.

---

## Contributing

Issues and pull requests are welcome. The architecture is intentionally
modular. Adding a backend, a new metric, or a new input format does not
require touching the core IR.

Run the test suite:

```bash
cargo test --workspace
cargo clippy --workspace
```

## Sponsoring

CleitonForge is developed and maintained by a solo independent researcher.
If this tool saved you time or helped you catch a bug, consider sponsoring:

**➜ [github.com/sponsors/cleitonaugusto](https://github.com/sponsors/cleitonaugusto)**

---

## License

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE).
