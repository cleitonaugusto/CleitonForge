# CleitonForge

**[English]** | [Português](README.pt-BR.md)

**CleitonForge** is an open-source, neutral benchmarking and interoperability
layer for quantum computing simulators, written in Rust. The CLI binary is
`cforge`.

---

## The problem

The quantum computing software ecosystem is fragmented across incompatible
frameworks — Qiskit (IBM), Cirq (Google), Amazon Braket, PennyLane (Xanadu)
— and several isolated Rust-based simulators (QuantRS2, qoqo, q1tsim, qvnt).
Each framework measures performance differently, uses different circuit formats,
and defines the same algorithms in incompatible ways.

The academic literature confirms the gap:

- A benchmark named "grover" can implement Grover's algorithm with a completely
  different oracle in each tool, making results incomparable across frameworks.
- Quantum tooling remains hard to integrate into CI/CD pipelines due to the
  absence of standardized interchange formats.
- Platform-agnostic benchmarking efforts exist (QED-C, MQT Bench) but their
  own authors acknowledge they are "a step", not a definitive solution.

## Why Rust, why now

- **OpenQASM 3** is the de facto interchange standard for the industry, but
  its reference implementation and most tooling are in Python.
- Rust-native OpenQASM parsers (`oq3_semantics`, `oq3_parser`, `oq3_syntax`)
  exist with 680k+ downloads each but have been unmaintained for over a year.
- Rust-based simulators exist but are isolated silos with no common API.
- Rust delivers predictable, low-overhead performance — critical for
  benchmarking where measurement noise must be minimized.

## Positioning

CleitonForge does not compete with Qiskit, Cirq, IBM, or any hardware vendor.
The goal is to be the **neutral layer** that all of them can plug into — an
impartial judge that measures and compares without favoring any backend.

- Never optimized to favor a specific simulator
- Tracks provenance of every gate (which framework, which original name)
- Plugin architecture (`SimulationBackend` trait) so any backend can be added
  by the community without touching the core

---

## Architecture

```
CleitonForge/                    Rust workspace
├── cforge-core/                 Canonical IR: Circuit, GateKind, Operation
├── cforge-parser/               OpenQASM 2 + OpenQASM 3 parsers
├── cforge-backends/             SimulationBackend trait + implementations
├── cforge-metrics/              Fidelity, depth, timing, memory measurement
├── cforge-cli/                  `cforge` binary (clap + comfy-table)
│   └── examples/
│       └── compare_grover.rs   Grover algorithm — pure Rust API example
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
- No external simulators required — all dependencies are pure Rust crates

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

**Planned:** additional backends (qoqo, q1tsim), noise modeling, Python
bindings (PyO3), extended OpenQASM 3 gate coverage, web dashboard.

---

## Contributing

Issues and pull requests are welcome. The architecture is intentionally
modular — adding a backend, a new metric, or a new input format does not
require touching the core IR.

Run the test suite:

```bash
cargo test --workspace
cargo clippy --workspace
```

---

## License

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE).
