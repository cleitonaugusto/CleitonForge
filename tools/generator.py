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

# Angles a compiler special-cases, plus the near-misses on either side of the
# thresholds it uses to decide.
#
# The default uniform draw over (-pi, pi) never lands on any of these: hitting an
# exact multiple of pi/2 has probability zero. That is not a small gap. #16594
# lives exactly at negative multiples of pi/2, and the only reason our campaigns
# ever found it is that sxdg hands over -pi/2 for free. We found it by accident,
# not by design, and every rotation bug that has no convenient gate handing it to
# us is invisible to a uniform generator.
#
# The thresholds are read from the source rather than guessed:
#   _CUTOFF_PRECISION = 1e-5 in commutation_cancellation.rs decides whether an
#   angle counts as a multiple of pi, so an angle 1e-6 away from pi/2 is treated
#   as pi/2 and synthesised into an exact gate.
#   ANGLE_ZERO_EPSILON = 1e-12 in euler_one_qubit_decomposer.rs decides whether a
#   rotation is dropped as zero.
#
# Sitting just inside and just outside each threshold is where a tolerance either
# rounds something it should not, or fails to round something it should.
_CUTOFF_PRECISION = 1e-5  # commutation_cancellation.rs:63
_ANGLE_ZERO_EPSILON = 1e-12  # euler_one_qubit_decomposer.rs:41

_SPECIAL_BASE = [
    0.0,
    math.pi / 8, math.pi / 4, math.pi / 2, 3 * math.pi / 4,
    math.pi, 3 * math.pi / 2, 2 * math.pi, 4 * math.pi,
]

_NEAR_MISS_DELTAS = [
    _CUTOFF_PRECISION * 0.1,   # inside the cutoff: gets rounded to the exact value
    _CUTOFF_PRECISION,         # on the cutoff itself
    _CUTOFF_PRECISION * 10,    # outside: must not be rounded
    _ANGLE_ZERO_EPSILON,       # below every threshold
]

# Exact special values, and near-misses, are kept apart on purpose. They answer
# different questions, and mixing them buries the answer to the first one.
#
# Exact values probe for real faults: rx(-3*pi/2) followed by sx accumulates to
# -pi and #16594 discards it, an error of 1.0 with no sxdg anywhere. That is a
# correctness bug, and a uniform generator cannot reach it.
#
# Near-misses probe whether the pass honours its own tolerance. An angle 1e-6
# from 3*pi/2 sits inside _CUTOFF_PRECISION, so the pass rounds it and emits an
# exact gate, introducing about 1e-6 of error. QCEC calls that not_equivalent,
# and it is doing its job: the difference is real. But it is the tolerance
# working as designed, not a fault, and it is the same class as the cp
# small-angle finding we catalogued as contract rather than bug.
#
# Measured on one run: mixing them gave 13 findings, of which 1 was the real bug
# and 12 were the tolerance. Turn near-misses on to ask about the tolerance, and
# leave them off when hunting faults.
SPECIAL_ANGLES = sorted(
    {sign * base for base in _SPECIAL_BASE for sign in (1.0, -1.0)}
)

NEAR_MISS_ANGLES = sorted(
    {sign * base + d
     for base in _SPECIAL_BASE
     for sign in (1.0, -1.0)
     for delta in _NEAR_MISS_DELTAS
     for d in (delta, -delta)}
)

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


def random_angle(rng: random.Random, boundary_prob: float = 0.5,
                 near_miss_prob: float = 0.0) -> float:
    """An angle, biased towards the values a compiler special-cases.

    Draws a near-miss with near_miss_prob, an exact special value with
    boundary_prob, and otherwise uniform over (-pi, pi). Keeping the uniform
    share matters: it is what exercises the generic path, and a generator that
    only ever produced special values would stop testing the ordinary one.
    """
    roll = rng.random()
    if roll < near_miss_prob:
        return rng.choice(NEAR_MISS_ANGLES)
    if roll < near_miss_prob + boundary_prob:
        return rng.choice(SPECIAL_ANGLES)
    return rng.uniform(-math.pi, math.pi)


def random_ops(rng: random.Random, num_qubits: int, depth: int,
               exclude: frozenset[str] = frozenset(),
               boundary_prob: float = 0.0,
               near_miss_prob: float = 0.0) -> list:
    """Random gate list [(name, qubits, params), ...].

    exclude drops gates from the pool. The use for it is suppressing a bug you
    have already reported so the search stops spending its budget rediscovering
    it: against Qiskit 2.5.0, excluding sxdg makes #16594 essentially
    unreachable, because that fault needs an accumulated X rotation at a
    negative multiple of pi/2, sx and x contribute positively, and random rx
    angles hit an exact multiple with probability zero.

    boundary_prob biases angles towards SPECIAL_ANGLES. It defaults to 0 so that
    existing campaigns keep their exact distribution and stay comparable; new
    campaigns should turn it on.
    """
    eligible = [g for g in GATES if g[1] <= num_qubits and g[0] not in exclude]
    if not eligible:
        raise ValueError("gate pool is empty after exclusions")
    weights = [g[3] for g in eligible]
    ops = []
    for _ in range(depth):
        name, nq, npar, _ = rng.choices(eligible, weights=weights)[0]
        qubits = rng.sample(range(num_qubits), nq)
        params = [random_angle(rng, boundary_prob, near_miss_prob) for _ in range(npar)]
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
    """QASM 2 for a gate list, with angles that round-trip.

    Angles print at full float precision (`repr`), not fixed decimals. A witness
    that does not reproduce is worthless, and rounding silently destroys the
    cases this generator exists to produce: `%.12f` turns 3*pi/2 + 1e-6 into
    something that reads as 3*pi/2, so the witness looks like an exact special
    angle and stops reproducing. That cost real debugging time once.
    """
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{num_qubits}];"]
    for name, qubits, params in ops:
        p = f"({','.join(repr(float(x)) for x in params)})" if params else ""
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
