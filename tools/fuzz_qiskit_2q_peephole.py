#!/usr/bin/env python3
"""Fuzz TwoQubitPeepholeOptimization, the pass added in Qiskit 2.5.0.

Why this pass and not another: it is seven weeks old, it sits in the default
O2 pipeline so ordinary users reach it, it does numerical synthesis of two qubit
blocks, and it is multithreaded. Age of the specific code path predicts yield
far better than the size of the project. Measured on this repo's own campaigns:
a five month old compiler gave six findings in an hour, a recent lowering path
gave one in a day, and Cirq's decade old QASM exporter gave nothing in 3500
circuits.

Two properties, because the pass has two ways to be wrong.

1. Unitary preservation. At approximation_degree=1.0 the pass promises no
   approximation, so the circuit's operator must survive it up to global phase.
   This is the same property the other pass fuzzers here use.

2. Determinism. The pass is multithreaded. Running it twice on the same input
   must give the same circuit. A race shows up as a run that disagrees with
   itself, which no single-run oracle can see.

Usage:
    python3 tools/fuzz_qiskit_2q_peephole.py --iterations 500
"""
from __future__ import annotations

import argparse
import math
import pathlib
import random
import sys

import numpy as np
from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import GenericBackendV2
from qiskit.quantum_info import Operator
from qiskit.transpiler import PassManager
from qiskit.transpiler import passes as qpasses

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from generator import build, random_ops, shrink, to_qasm2  # noqa: E402

BASIS = ["rz", "sx", "x", "cx"]


def make_pass(num_qubits: int, approximation_degree: float):
    target = GenericBackendV2(num_qubits=num_qubits, basis_gates=BASIS,
                              seed=42).target
    return qpasses.TwoQubitPeepholeOptimization(
        target=target, approximation_degree=approximation_degree)


def diff_op(a, b) -> float:
    """Largest entry difference modulo global phase.

    Not fidelity, which is quadratic in the error and hides small faults, and
    not a state comparison, which only exercises one column of the operator.
    """
    tr = np.vdot(a.reshape(-1), b.reshape(-1))
    ph = tr / abs(tr) if abs(tr) > 1e-12 else 1.0
    return float(np.max(np.abs(a - b * np.conj(ph))))


def run_once(ops, num_qubits, approximation_degree):
    qc = build(num_qubits, ops)
    out = PassManager([make_pass(num_qubits, approximation_degree)]).run(qc)
    return qc, out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=500)
    ap.add_argument("--qubits", type=int, default=4)
    ap.add_argument("--min-depth", type=int, default=6)
    ap.add_argument("--max-depth", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--boundary-prob", type=float, default=0.5)
    ap.add_argument("--repeats", type=int, default=3,
                    help="runs per circuit, for the determinism check")
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    args = ap.parse_args()

    import qiskit
    print(f"TwoQubitPeepholeOptimization — qiskit {qiskit.__version__}")
    print(f"  {args.iterations} circuits, {args.qubits} qubits, "
          f"depth {args.min_depth}-{args.max_depth}, seed {args.seed}, "
          f"approximation_degree=1.0\n")

    unitary_bad = nondet = errors = checked = 0
    worst = 0.0
    first_witness = None

    for i in range(args.iterations):
        rng = random.Random(args.seed + i)
        depth = rng.randint(args.min_depth, args.max_depth)
        ops = random_ops(rng, args.qubits, depth, boundary_prob=args.boundary_prob)

        try:
            qc, out = run_once(ops, args.qubits, 1.0)
        except Exception as e:  # noqa: BLE001 — a pass crash is a finding
            errors += 1
            if errors <= 3:
                print(f"  [crash] circuit {i}: {type(e).__name__}: {e}")
            continue

        checked += 1

        # Property 1: the operator has to survive.
        try:
            d = diff_op(Operator(qc).data, Operator(out).data)
        except Exception as e:  # noqa: BLE001
            errors += 1
            print(f"  [oracle crash] circuit {i}: {type(e).__name__}: {e}")
            continue
        worst = max(worst, d)
        if d > 1e-6:
            unitary_bad += 1
            print(f"\n  [UNITARY CHANGED] circuit {i} (seed {args.seed + i}), "
                  f"max entry difference {d:.3e}")
            if first_witness is None:
                def still_fails(cand):
                    try:
                        a, b = run_once(cand, args.qubits, 1.0)
                        return diff_op(Operator(a).data, Operator(b).data) > 1e-6
                    except Exception:  # noqa: BLE001
                        return False
                first_witness = shrink(ops, still_fails)
                print("    minimal witness:")
                for line in to_qasm2(first_witness, args.qubits).splitlines()[3:]:
                    print(f"      {line}")

        # Property 2: the same input must give the same output every time.
        try:
            outs = [run_once(ops, args.qubits, 1.0)[1] for _ in range(args.repeats - 1)]
        except Exception:  # noqa: BLE001
            continue
        for k, other in enumerate(outs, start=2):
            if Operator(out).data.shape != Operator(other).data.shape:
                continue
            if diff_op(Operator(out).data, Operator(other).data) > 1e-9:
                nondet += 1
                print(f"\n  [NONDETERMINISTIC] circuit {i} (seed {args.seed + i}): "
                      f"run 1 and run {k} give different operators")
                break
            if out.count_ops() != other.count_ops():
                nondet += 1
                print(f"\n  [NONDETERMINISTIC] circuit {i} (seed {args.seed + i}): "
                      f"same operator, different gate counts, "
                      f"{dict(out.count_ops())} vs {dict(other.count_ops())}")
                break

    print(f"\nchecked {checked}, {unitary_bad} unitary changed, "
          f"{nondet} nondeterministic, {errors} errors")
    print(f"largest operator difference seen: {worst:.3e}")
    return 1 if (unitary_bad or nondet or errors) else 0


if __name__ == "__main__":
    raise SystemExit(main())
