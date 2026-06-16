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

## Status

Early development. The canonical IR, parsers, backends, metrics, and CLI
are being built incrementally — see the crate list above for current
scope.

## License

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE).
