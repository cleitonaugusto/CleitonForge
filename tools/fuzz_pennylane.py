#!/usr/bin/env python3
"""Phase-sensitive differential fuzzing of PennyLane — Phase 3.

Mirrors tools/fuzz_qiskit.py. Two devices simulate the same circuit and
their statevectors are compared modulo global phase:

    lightning.qubit (C++ Kokkos backend)  vs  default.qubit (Python)

lightning.qubit is the performance backend most PennyLane users run in
production; default.qubit is the reference implementation. Any N1
divergence between them is a backend bug (PennyLane issue #9837 tracks
this class). Witnesses shrink to a 1-minimal core and go to the bug zoo.

Usage:
    python3 tools/fuzz_pennylane.py --iterations 3000 --seed 42
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

TOL = 1e-9

# (PennyLane op, n_qubits, n_params, weight) — matched to fuzz_qiskit.
# Each entry is (name, builder) where builder(params, wires) applies it.
GATES = [
    ("RZ", 1, 1, 5), ("PhaseShift", 1, 1, 5), ("CRZ", 2, 1, 5), ("ControlledPhaseShift", 2, 1, 5),
    ("S", 1, 0, 5), ("Sdg", 1, 0, 5), ("T", 1, 0, 5), ("Tdg", 1, 0, 5),
    ("SX", 1, 0, 5), ("SXdg", 1, 0, 5),
    ("Hadamard", 1, 0, 4),
    ("RX", 1, 1, 2), ("RY", 1, 1, 2),
    ("CNOT", 2, 0, 3), ("CY", 2, 0, 1), ("CZ", 2, 0, 1), ("CH", 2, 0, 1),
    ("CRX", 2, 1, 1), ("CRY", 2, 1, 1),
    ("PauliX", 1, 0, 1), ("PauliY", 1, 0, 1), ("PauliZ", 1, 0, 1),
    ("SWAP", 2, 0, 1), ("Toffoli", 3, 0, 1), ("CSWAP", 3, 0, 1),
]


def apply(name, params, wires):
    if name == "Sdg":
        qml.adjoint(qml.S)(wires=wires[0])
    elif name == "Tdg":
        qml.adjoint(qml.T)(wires=wires[0])
    elif name == "SXdg":
        qml.adjoint(qml.SX)(wires=wires[0])
    elif name in ("S", "T", "SX", "Hadamard", "PauliX", "PauliY", "PauliZ"):
        getattr(qml, name)(wires=wires[0])
    elif name in ("RZ", "RX", "RY", "PhaseShift"):
        getattr(qml, name)(params[0], wires=wires[0])
    elif name in ("CNOT", "CY", "CZ", "CH", "SWAP"):
        getattr(qml, name)(wires=wires)
    elif name in ("CRZ", "CRX", "CRY", "ControlledPhaseShift"):
        getattr(qml, name)(params[0], wires=wires)
    elif name in ("Toffoli", "CSWAP"):
        getattr(qml, name)(wires=wires)
    else:
        raise ValueError(f"unhandled gate {name}")


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


def statevector(ops, num_qubits, device):
    dev = qml.device(device, wires=num_qubits)

    @qml.qnode(dev)
    def circuit():
        for name, wires, params in ops:
            apply(name, params, wires)
        return qml.state()

    return np.asarray(circuit())


def n1_amplitude(a, b):
    overlap = np.vdot(a, b)
    phase = overlap / abs(overlap) if abs(overlap) > 1e-15 else 1.0
    return float(np.max(np.abs(a * phase - b)))


def divergence(ops, num_qubits):
    try:
        a = statevector(ops, num_qubits, "default.qubit")
        b = statevector(ops, num_qubits, "lightning.qubit")
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
            candidate = current[:i] + current[i + 1:]
            if candidate and still_fails(candidate):
                current = candidate
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
    ap.add_argument("--iterations", type=int, default=3000)
    ap.add_argument("--max-qubits", type=int, default=4)
    ap.add_argument("--max-depth", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    args = ap.parse_args()

    print(f"fuzz_pennylane — pennylane {qml.__version__}")
    print(f"  differential: lightning.qubit vs default.qubit")
    print(f"  iterations: {args.iterations}  seed: {args.seed}  tol: {args.tol:g}\n")

    found = 0
    worst = 0.0
    for i in range(args.iterations):
        rng = random.Random(args.seed + i)
        nq = rng.randint(1, args.max_qubits)
        depth = rng.randint(1, args.max_depth)
        ops = random_ops(rng, nq, depth)

        dist = divergence(ops, nq)
        if dist is None:
            continue
        worst = max(worst, dist)
        if dist <= args.tol:
            continue

        found += 1
        minimal = shrink(ops, nq, lambda cand: (divergence(cand, nq) or 0.0) > args.tol)
        final = divergence(minimal, nq) or 0.0
        print(f"── divergence #{found} [lightning-vs-default] circuit {i}, N1 = {dist:.3e}")
        print(f"  minimal witness ({len(minimal)} gates):")
        for line in describe(minimal):
            print(f"    {line}")

        entry = {
            "id": f"pennylane-lightning-{found:03}",
            "target": "PennyLane lightning.qubit vs default.qubit",
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

        if (i + 1) % 500 == 0:
            print(f"[{i + 1}/{args.iterations}] worst N1: {worst:.2e}")

    print(f"\nResult: {found} divergence(s) in {args.iterations} circuits.")
    print(f"Worst N1 distance: {worst:.3e}  ({'CLEAN' if worst <= args.tol else 'BUG'})")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
