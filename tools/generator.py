#!/usr/bin/env python3
"""Weighted random circuit generator.

Extracted from fuzz_qiskit.py, where it had been copied into six harnesses.
This is the part of the fuzzer that is actually ours: QCEC is a better oracle
than anything we would write, but nobody else weights a gate pool for
phase-convention bugs.

The weights are the whole point. Phase-family gates (rz, p, crz, cp) and the
discrete complex gates (s, sdg, t, tdg, sx, sxdg) get 5x, so that circuits
mixing a possibly-bugged rotation with a known-good complex gate dominate the
search. A uniform pool wastes most of its budget on Cliffords, where sign
errors cancel.

sxdg carries its own history: it is the gate that exposed Qiskit #16594.
"""

from __future__ import annotations

import math
import random

# (name, n_qubits, n_params, weight)
GATES = [
    ("rz", 1, 1, 5), ("p", 1, 1, 5), ("crz", 2, 1, 5), ("cp", 2, 1, 5),
    ("s", 1, 0, 5), ("sdg", 1, 0, 5), ("t", 1, 0, 5), ("tdg", 1, 0, 5),
    ("sx", 1, 0, 5), ("sxdg", 1, 0, 5),
    ("h", 1, 0, 4),
    ("rx", 1, 1, 2), ("ry", 1, 1, 2), ("u", 1, 3, 2),
    ("cx", 2, 0, 3), ("cy", 2, 0, 1), ("cz", 2, 0, 1), ("ch", 2, 0, 1),
    ("csx", 2, 0, 2), ("crx", 2, 1, 1), ("cry", 2, 1, 1),
    ("x", 1, 0, 1), ("y", 1, 0, 1), ("z", 1, 0, 1),
    ("swap", 2, 0, 1), ("ccx", 3, 0, 1), ("cswap", 3, 0, 1),
]


def random_ops(rng: random.Random, num_qubits: int, depth: int) -> list:
    """Random gate list [(name, qubits, params), ...]."""
    eligible = [g for g in GATES if g[1] <= num_qubits]
    weights = [g[3] for g in eligible]
    ops = []
    for _ in range(depth):
        name, nq, npar, _ = rng.choices(eligible, weights=weights)[0]
        qubits = rng.sample(range(num_qubits), nq)
        params = [rng.uniform(-math.pi, math.pi) for _ in range(npar)]
        ops.append((name, qubits, params))
    return ops


def build(num_qubits: int, ops):
    """Materialize a gate list into a QuantumCircuit."""
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(num_qubits)
    for name, qubits, params in ops:
        getattr(qc, name)(*params, *qubits)
    return qc


def to_qasm2(ops, num_qubits: int) -> str:
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{num_qubits}];"]
    for name, qubits, params in ops:
        p = f"({','.join(f'{x:.12f}' for x in params)})" if params else ""
        q = ",".join(f"q[{i}]" for i in qubits)
        lines.append(f"{name}{p} {q};")
    return "\n".join(lines) + "\n"


def shrink(ops, still_fails) -> list:
    """Greedy 1-minimal reduction by gate removal.

    still_fails(candidate_ops) -> bool. Oracle-agnostic on purpose: it works
    with the QCEC verdict just as it did with the amplitude distance.
    """
    current = list(ops)
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(current):
            candidate = current[:i] + current[i + 1:]
            if candidate and still_fails(candidate):
                current = candidate
                changed = True
            else:
                i += 1
    return current
