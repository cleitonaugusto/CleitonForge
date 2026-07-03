# CleitonForge

**CleitonForge** is an open-source benchmarking and interoperability layer
for quantum computing simulators, written in Rust. The CLI binary is
`cforge`.

## The problem

The quantum computing ecosystem is fragmented across incompatible
frameworks — Qiskit (IBM), Cirq (Google), Amazon Braket, PennyLane
(Xanadu) — and several isolated Rust-based simulators (QuantRS2, qoqo,
q1tsim). Each has its own way of measuring performance and correctness.
The academic literature confirms the gap:

- There is no consolidated agreement on which benchmarks to use for
  empirical evaluation of quantum tools — a benchmark called "grover" can
  implement Grover's algorithm with a different oracle in every tool,
  making results incomparable across frameworks.
- Quantum tooling remains hard to integrate into IDEs and CI/CD pipelines
  due to a lack of standardized formats and real interoperability.
- Platform-agnostic benchmarking efforts exist (QED-C, MQT Bench), but
  their own authors acknowledge they don't fully solve the problem — they
  are "a step", not a definitive solution.

## Why Rust, why now

- OpenQASM 3 is already the de facto interchange standard in the
  industry, but its reference implementation is in Python.
- OpenQASM parsers already exist in Rust (`oq3_parser`, `oq3_lexer`,
  `oq3_semantics`, `oq3_syntax`) with over 680k downloads each, but have
  not been updated in over a year — a sign of opportunity to extend, not
  a reason to start from scratch.
- Rust-based simulators exist but are fragmented and isolated: QuantRS2,
  qoqo, q1tsim, qvnt — none has a unified benchmarking layer.
- Rust offers superior performance over Python for heavy structural
  analysis of large circuits.

## Positioning

CleitonForge does not compete with Qiskit, Cirq, IBM, or any hardware
vendor. The goal is to be the **neutral layer** that all of them can use —
an impartial judge that measures and compares, without favoring any
backend. This means:

- Never optimizing the tool to favor a specific simulator
- Always crediting and tracking the provenance of each source framework
- A plugin/trait architecture so any new backend can be added by the
  community without rewriting the core

## Architecture

```
cleitonforge/                    (Rust workspace)
├── cforge-core/                 (canonical IR + shared types)
├── cforge-parser/               (multi-format input layer: OpenQASM 2/3)
├── cforge-backends/             (pluggable SimulationBackend trait + impls)
├── cforge-metrics/              (standardized metrics computation)
├── cforge-cli/                  (command-line interface, binary: cforge)
└── examples/                    (end-to-end usage examples)
```

## Quick Start

### Prerequisites

- Rust 1.96+ (`rustup update stable`)
- No external simulators required — everything runs inside the workspace

### Build

```bash
git clone https://github.com/cleitonaugusto/CleitonForge.git
cd CleitonForge
cargo build --release
```

The binary lands at `target/release/cforge`.

### Run a circuit on both backends

```bash
# Bell state — compare statevector-native vs statevector-quantrs2
cargo run --release --bin cforge -- run \
    --circuit examples/bell.qasm \
    --backends statevector,quantrs2 \
    --shots 1024
```

Expected output:

```
Circuit: 2 qubits  |  2 gates  |  depth 2
┌──────────────────────┬───────────┬───────┬───────┬──────────┬───────┐
│ Backend              │ Time (ms) │ Depth │ Gates │ Fidelity │ Shots │
╞══════════════════════╪═══════════╪═══════╪═══════╪══════════╪═══════╡
│ statevector-native   │ 0.002     │ 2     │ 2     │ 1.000000 │ 1024  │
├──────────────────────┼───────────┼───────┼───────┼──────────┼───────┤
│ statevector-quantrs2 │ 0.004     │ 2     │ 2     │ 1.000000 │ 1024  │
└──────────────────────┴───────────┴───────┴───────┴──────────┴───────┘
```

### Validate a circuit (no simulation)

```bash
cargo run --bin cforge -- validate --circuit examples/bell.qasm
```

### Grover search end-to-end example

This example builds a 3-qubit Grover circuit entirely in Rust (no `.qasm`
file), searches for the state |101⟩, runs it on both backends, and
confirms they agree:

```bash
cargo run --example compare_grover -p cforge-cli
```

Expected output (abridged):

```
CleitonForge — Grover search example
Target state : |101⟩  (index 5, q0=1 q1=0 q2=1)
Qubits       : 3   (N=8 states)
Iterations   : 2   (sin²(5θ) ≈ 94.8 % expected)

Circuit      : 43 gates  |  depth 21

Backend  : statevector-native
  Top state  : |101⟩  index 5  prob = 0.9453 (94.5 %)
  Fidelity   : 1.00000000
  Top counts (1024 shots):
    |101⟩    965 shots  (94.2 %)

Backend  : statevector-quantrs2
  Top state  : |101⟩  index 5  prob = 0.9453 (94.5 %)
  Fidelity   : 1.00000000
  Top counts (1024 shots):
    |101⟩    965 shots  (94.2 %)

Cross-backend fidelity (native vs quantrs2): 1.00000000
Both backends agree: YES ✓
```

The measured 94.5 % probability matches the theoretical prediction of
sin²(5θ) ≈ 94.8 % for N=8, M=1, k=2 Grover iterations — confirming
correctness of both the circuit and the two independent backends.

### Supported backends

| Name          | Flag           | Description                                      |
|---------------|----------------|--------------------------------------------------|
| Native SV     | `statevector`  | Custom state-vector sim built inside CleitonForge |
| QuantRS2      | `quantrs2`     | Uses `quantrs2-core` gate matrices               |

### Supported input formats

| Format       | Notes                              |
|--------------|------------------------------------|
| OpenQASM 2.0 | Auto-detected; no stdgates include needed |
| OpenQASM 3   | Auto-detected via `OPENQASM 3` header |

## Status

All core phases are complete and the CLI is functional:

| Phase | Crate              | Status   |
|-------|--------------------|----------|
| 0     | workspace setup    | ✓ done   |
| 1     | `cforge-core`      | ✓ done   |
| 2     | `cforge-parser`    | ✓ done   |
| 3     | `cforge-backends`  | ✓ done   |
| 4     | `cforge-metrics`   | ✓ done   |
| 5     | `cforge-cli`       | ✓ done   |
| 6     | examples + docs    | ✓ done   |

Planned next: additional backends (qoqo, q1tsim), noise modeling,
OpenQASM 3 extended gate coverage, web UI.

## License

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE).
