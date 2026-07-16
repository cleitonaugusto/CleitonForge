#!/usr/bin/env python3
"""Metamorphic fuzzing of PennyLane's optimization transforms — Phase 3.

The Qiskit bug we found lived in a *cancellation pass*, not in the raw
simulator. PennyLane ships the same class of passes as circuit
transforms:

    cancel_inverses, merge_rotations, single_qubit_fusion,
    commute_controlled, pattern_matching_optimization, undo_swaps

Each must preserve the statevector modulo global phase. This harness
builds a random circuit, applies one transform, and compares the state
before and after with the N1 amplitude oracle. Any divergence is a
soundness bug in that transform. Witnesses shrink to a 1-minimal core.

Usage:
    python3 tools/fuzz_pennylane_transforms.py --iterations 4000 --seed 42
    python3 tools/fuzz_pennylane_transforms.py --transform cancel_inverses
"""

import argparse
import json
import pathlib
import random
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

import pennylane as qml  # noqa: E402
from pennylane.tape import QuantumTape  # noqa: E402

TOL = 1e-9

TRANSFORMS = {
    "cancel_inverses": qml.transforms.cancel_inverses,
    "merge_rotations": qml.transforms.merge_rotations,
    "single_qubit_fusion": qml.transforms.single_qubit_fusion,
    "commute_controlled": qml.transforms.commute_controlled,
    "undo_swaps": qml.transforms.undo_swaps,
}

# Matched to fuzz_pennylane.py; sxdg/sx-heavy to probe the analogue of the
# Qiskit CommutativeCancellation bug.
GATES = [
    ("RZ", 1, 1, 5), ("PhaseShift", 1, 1, 5), ("CRZ", 2, 1, 5), ("ControlledPhaseShift", 2, 1, 5),
    ("S", 1, 0, 5), ("Sdg", 1, 0, 5), ("T", 1, 0, 5), ("Tdg", 1, 0, 5),
    ("SX", 1, 0, 5), ("SXdg", 1, 0, 5),
    ("Hadamard", 1, 0, 4),
    ("RX", 1, 1, 2), ("RY", 1, 1, 2),
    ("CNOT", 2, 0, 3), ("CY", 2, 0, 1), ("CZ", 2, 0, 1),
    ("CRX", 2, 1, 1), ("CRY", 2, 1, 1),
    ("PauliX", 1, 0, 1), ("PauliY", 1, 0, 1), ("PauliZ", 1, 0, 1),
    ("SWAP", 2, 0, 1), ("Toffoli", 3, 0, 1),
]


def make_op(name, params, wires):
    if name == "Sdg":
        return qml.adjoint(qml.S(wires=wires[0]))
    if name == "Tdg":
        return qml.adjoint(qml.T(wires=wires[0]))
    if name == "SXdg":
        return qml.adjoint(qml.SX(wires=wires[0]))
    if name in ("S", "T", "SX", "Hadamard", "PauliX", "PauliY", "PauliZ"):
        return getattr(qml, name)(wires=wires[0])
    if name in ("RZ", "RX", "RY", "PhaseShift"):
        return getattr(qml, name)(params[0], wires=wires[0])
    if name in ("CNOT", "CY", "CZ", "SWAP"):
        return getattr(qml, name)(wires=wires)
    if name in ("CRZ", "CRX", "CRY", "ControlledPhaseShift"):
        return getattr(qml, name)(params[0], wires=wires)
    if name == "Toffoli":
        return qml.Toffoli(wires=wires)
    raise ValueError(name)


def random_ops(rng, num_qubits, depth):
    eligible = [g for g in GATES if g[1] <= num_qubits]
    weights = [g[3] for g in eligible]
    ops = []
    for _ in range(depth):
        name, nq, npar, _ = rng.choices(eligible, weights=weights)[0]
        wires = rng.sample(range(num_qubits), nq)
        params = [rng.uniform(-np.pi, np.pi) for _ in range(npar)]
        ops.append((name, wires, params))
    return ops


def state_of(ops, num_qubits):
    with QuantumTape() as tape:
        for name, wires, params in ops:
            make_op(name, params, wires)
    dev = qml.device("default.qubit", wires=num_qubits)
    program = dev.preprocess_transforms()
    tapes, _ = program([tape])
    res = dev.execute(tapes)
    # statevector via qml.state on the same circuit
    return _statevector(ops, num_qubits)


def _statevector(ops, num_qubits):
    dev = qml.device("default.qubit", wires=num_qubits)

    @qml.qnode(dev)
    def circ():
        for name, wires, params in ops:
            make_op(name, params, wires)
        return qml.state()

    return np.asarray(circ())


def _statevector_transformed(ops, num_qubits, transform):
    with QuantumTape() as tape:
        for name, wires, params in ops:
            make_op(name, params, wires)
    new_tapes, _ = transform(tape)
    tape2 = new_tapes[0]

    dev = qml.device("default.qubit", wires=num_qubits)

    @qml.qnode(dev)
    def circ():
        for op in tape2.operations:
            qml.apply(op)
        return qml.state()

    return np.asarray(circ())


def n1_amplitude(a, b):
    overlap = np.vdot(a, b)
    phase = overlap / abs(overlap) if abs(overlap) > 1e-15 else 1.0
    return float(np.max(np.abs(a * phase - b)))


def divergence(ops, num_qubits, transform):
    try:
        a = _statevector(ops, num_qubits)
        b = _statevector_transformed(ops, num_qubits, transform)
        if a.shape != b.shape:
            return float("inf")
        return n1_amplitude(a, b)
    except Exception as e:  # noqa: BLE001
        print(f"  [crash] {type(e).__name__}: {e}", file=sys.stderr)
        return None


def shrink(ops, num_qubits, still_fails):
    current = list(ops)
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(current):
            cand = current[:i] + current[i + 1:]
            if cand and still_fails(cand):
                current = cand
                changed = True
            else:
                i += 1
    return current


def describe(ops):
    out = []
    for name, wires, params in ops:
        p = f"({','.join(f'{x:.4f}' for x in params)})" if params else ""
        out.append(f"{name}{p} {wires}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=4000)
    ap.add_argument("--max-qubits", type=int, default=4)
    ap.add_argument("--max-depth", type=int, default=12)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--transform", choices=list(TRANSFORMS) + ["all"], default="all")
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    args = ap.parse_args()

    which = list(TRANSFORMS) if args.transform == "all" else [args.transform]
    print(f"fuzz_pennylane_transforms — pennylane {qml.__version__}")
    print(f"  transforms: {which}")
    print(f"  iterations: {args.iterations}  seed: {args.seed}  tol: {args.tol:g}\n")

    total_found = 0
    for tname in which:
        transform = TRANSFORMS[tname]
        found = 0
        worst = 0.0
        for i in range(args.iterations):
            rng = random.Random(args.seed + i)
            nq = rng.randint(1, args.max_qubits)
            depth = rng.randint(1, args.max_depth)
            ops = random_ops(rng, nq, depth)

            d = divergence(ops, nq, transform)
            if d is None:
                continue
            worst = max(worst, d)
            if d <= args.tol:
                continue

            found += 1
            total_found += 1
            minimal = shrink(ops, nq, lambda c: (divergence(c, nq, transform) or 0.0) > args.tol)
            final = divergence(minimal, nq, transform) or 0.0
            print(f"── [{tname}] divergence #{found} circuit {i}, N1 = {d:.3e}")
            print(f"  minimal witness ({len(minimal)} gates):")
            for line in describe(minimal):
                print(f"    {line}")

            entry = {
                "id": f"pennylane-{tname}-{found:03}",
                "target": f"PennyLane {tname} transform",
                "pennylane_version": qml.__version__,
                "generator_seed": args.seed + i,
                "num_qubits": nq,
                "gates": [{"gate": g, "wires": w, "params": p} for g, w, p in minimal],
                "n1_amplitude_distance": final,
            }
            args.zoo_dir.mkdir(exist_ok=True)
            out = args.zoo_dir / f"{entry['id']}.json"
            out.write_text(json.dumps(entry, indent=2))
            print(f"  zoo entry: {out}\n")

        tag = "CLEAN" if worst <= args.tol else "BUG"
        print(f"[{tname}] {found} divergence(s), worst N1 = {worst:.3e}  ({tag})\n")

    print(f"Total divergences across transforms: {total_found}")
    return 1 if total_found else 0


if __name__ == "__main__":
    sys.exit(main())
