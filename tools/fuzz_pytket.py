#!/usr/bin/env python3
"""Phase-sensitive metamorphic fuzzing of pytket (tket) — Phase 3.

Mirrors tools/fuzz_qiskit.py, but the target is Cambridge Quantum's
tket optimizer. The metamorphic property: a semantics-preserving
optimization pass must not change the statevector modulo global phase.

    |psi(C)>  vs  |psi(FullPeepholeOptimise(C))>

Any N1 divergence is a soundness bug in the optimizer. Witnesses are
shrunk to a 1-minimal core and written to the bug zoo.

Note on units: tket parameterizes rotations in half-turns, so Rz(a)
rotates by a*pi. The generator samples params in [-2, 2] half-turns.
This is internally consistent because both circuits (original and
optimized) live in tket, so the units never leave the framework.

Usage:
    python3 tools/fuzz_pytket.py --iterations 3000 --seed 42
"""

import argparse
import json
import pathlib
import random
import sys

import numpy as np
from pytket import Circuit
from pytket.passes import FullPeepholeOptimise

TOL = 1e-9

# (method name on Circuit, n_qubits, n_params, weight). Chosen to match
# the fuzz_qiskit pool: phase family and discrete complex gates 5x.
GATES = [
    ("Rz", 1, 1, 5), ("U1", 1, 1, 5), ("CRz", 2, 1, 5), ("CU1", 2, 1, 5),
    ("S", 1, 0, 5), ("Sdg", 1, 0, 5), ("T", 1, 0, 5), ("Tdg", 1, 0, 5),
    ("SX", 1, 0, 5), ("SXdg", 1, 0, 5),
    ("H", 1, 0, 4),
    ("Rx", 1, 1, 2), ("Ry", 1, 1, 2),
    ("CX", 2, 0, 3), ("CY", 2, 0, 1), ("CZ", 2, 0, 1), ("CH", 2, 0, 1),
    ("CSX", 2, 0, 2), ("CRx", 2, 1, 1), ("CRy", 2, 1, 1),
    ("X", 1, 0, 1), ("Y", 1, 0, 1), ("Z", 1, 0, 1),
    ("SWAP", 2, 0, 1), ("CCX", 3, 0, 1), ("CSWAP", 3, 0, 1),
]


def random_ops(rng, num_qubits, depth):
    eligible = [g for g in GATES if g[1] <= num_qubits]
    weights = [g[3] for g in eligible]
    ops = []
    for _ in range(depth):
        name, nq, npar, _ = rng.choices(eligible, weights=weights)[0]
        qubits = rng.sample(range(num_qubits), nq)
        params = [rng.uniform(-2.0, 2.0) for _ in range(npar)]
        ops.append((name, qubits, params))
    return ops


def build(num_qubits, ops):
    c = Circuit(num_qubits)
    for name, qubits, params in ops:
        getattr(c, name)(*params, *qubits)
    return c


def n1_amplitude(a, b):
    overlap = np.vdot(a, b)
    phase = overlap / abs(overlap) if abs(overlap) > 1e-15 else 1.0
    return float(np.max(np.abs(a * phase - b)))


def optimise_divergence(ops, num_qubits):
    """N1 distance between a circuit and its FullPeepholeOptimise form."""
    try:
        c = build(num_qubits, ops)
        sv = c.get_statevector()
        c2 = c.copy()
        FullPeepholeOptimise().apply(c2)
        sv2 = c2.get_statevector()
        return n1_amplitude(sv, sv2)
    except Exception as e:  # noqa: BLE001 — a crash is also a finding
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
    parts = []
    for name, qubits, params in ops:
        p = f"({','.join(f'{x:.4f}' for x in params)})" if params else ""
        parts.append(f"{name}{p} {qubits}")
    return parts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=3000)
    ap.add_argument("--max-qubits", type=int, default=4)
    ap.add_argument("--max-depth", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    args = ap.parse_args()

    import pytket

    print(f"fuzz_pytket — pytket {pytket.__version__}")
    print(f"  target: FullPeepholeOptimise soundness (metamorphic)")
    print(f"  iterations: {args.iterations}  seed: {args.seed}  tol: {args.tol:g}\n")

    found = 0
    worst = 0.0
    for i in range(args.iterations):
        rng = random.Random(args.seed + i)
        nq = rng.randint(1, args.max_qubits)
        depth = rng.randint(1, args.max_depth)
        ops = random_ops(rng, nq, depth)

        dist = optimise_divergence(ops, nq)
        if dist is None:
            continue
        worst = max(worst, dist)
        if dist <= args.tol:
            continue

        found += 1
        minimal = shrink(ops, nq, lambda cand: (optimise_divergence(cand, nq) or 0.0) > args.tol)
        final = optimise_divergence(minimal, nq) or 0.0
        print(f"── divergence #{found} [pytket-FullPeepholeOptimise] circuit {i}, N1 = {dist:.3e}")
        print(f"  minimal witness ({len(minimal)} gates):")
        for line in describe(minimal):
            print(f"    {line}")

        entry = {
            "id": f"pytket-fullpeephole-{found:03}",
            "target": "pytket (tket) FullPeepholeOptimise",
            "pytket_version": pytket.__version__,
            "generator_seed": args.seed + i,
            "num_qubits": nq,
            "gates": [{"gate": g, "qubits": qs, "params": ps} for g, qs, ps in minimal],
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
