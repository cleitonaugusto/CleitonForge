# CleitonForge

**Neutral benchmarking layer for quantum circuit simulators — Python bindings**

[![PyPI](https://img.shields.io/pypi/v/cleitonforge)](https://pypi.org/project/cleitonforge/)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue)](https://github.com/cleitonaugusto/CleitonForge/blob/main/LICENSE)
[![GitHub](https://img.shields.io/badge/GitHub-CleitonForge-black)](https://github.com/cleitonaugusto/CleitonForge)

**CleitonForge** is an open-source, Rust-powered neutral benchmarking and
interoperability layer for quantum computing simulators, created by
**Cleiton Augusto Correa Bezerra**.

This package (`cleitonforge`) exposes the full CleitonForge engine to Python —
build circuits, run simulations, parse OpenQASM, and measure performance metrics,
all from Python with the speed of native Rust.

---

## Install

```bash
pip install cleitonforge
```

No Rust toolchain required — pre-built wheels for Linux, macOS and Windows.

---

## Quick start

### Build and run a circuit

```python
import cleitonforge as cf

# Bell state
circuit = cf.Circuit(2)
circuit.h(0)
circuit.cx(0, 1)

result = cf.run(circuit, shots=1024)
print(result.counts)       # {'00': 512, '11': 512}
print(result.top_states(2))  # [('11', 0.5), ('00', 0.5)]
```

### Run from OpenQASM

```python
qasm = """
OPENQASM 2.0;
qreg q[2];
h q[0];
cx q[0],q[1];
"""

result = cf.run_qasm(qasm, shots=1024)
print(result.counts)
```

### Validate a circuit without simulating

```python
info = cf.validate_qasm(qasm)
# {'qubits': 2, 'gates': 2, 'depth': 2}
```

### Reproducible results with seed

```python
r1 = cf.run(circuit, shots=1024, seed=42)
r2 = cf.run(circuit, shots=1024, seed=42)
assert r1.counts == r2.counts  # always True
```

---

## Supported gates

Full OpenQASM 3 `stdgates.inc` gate set: `h`, `x`, `y`, `z`, `s`, `sdg`,
`t`, `tdg`, `sx`, `sxdg`, `rx`, `ry`, `rz`, `p`, `u`, `cx`, `cy`, `cz`,
`ch`, `csx`, `crx`, `cry`, `crz`, `cp`, `cu`, `swap`, `ccx`, `cswap`.

IBM/Qiskit aliases (`id`, `cnot`, `toffoli`, `fredkin`, `u2`, `u3`) are
automatically translated.

---

## Backends

| Name | Description |
|------|-------------|
| `"statevector"` | Native Rust state-vector simulator (default) |
| `"quantrs2"` | [QuantRS2](https://github.com/Qhronoboros/quantrs) backend |

```python
result = cf.run(circuit, shots=1024, backend="quantrs2")
```

---

## Statevector access

```python
result = cf.run(circuit, shots=0)  # shots=0 → exact statevector
for re, im in result.statevector:
    print(complex(re, im))
```

---

## Why CleitonForge

The quantum computing ecosystem is fragmented: Qiskit, Cirq, PennyLane, Amazon
Braket, and multiple Rust simulators each define benchmarks differently.
CleitonForge is the **neutral arbiter** — it measures and compares any backend
using the same canonical IR, without favoring any vendor.

- Canonical IR: every gate carries provenance (source framework, original name)
- Plugin architecture: add any backend via the `SimulationBackend` Rust trait
- OpenQASM 2 and 3 parsing natively in Rust
- Metrics over canonical IR guarantee fair cross-backend comparison

---

## CLI tool

The Rust CLI (`cforge`) is available in the main repository:

```bash
cforge run --circuit bell.qasm --backends statevector,quantrs2
cforge validate --circuit circuit.qasm
```

---

## Author

**Cleiton Augusto Correa Bezerra**
Brazilian systems analyst and quantum computing researcher.
[GitHub](https://github.com/cleitonaugusto) · [Email](mailto:augusto.cleiton@gmail.com)

---

## License

Apache 2.0 — same license as OpenQASM and the Rust ecosystem.
