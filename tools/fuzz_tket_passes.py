#!/usr/bin/env python3
"""Metamorphic fuzzing of tket optimization passes — Phase 3.

Same idea as fuzz_pennylane_transforms.py, for Cambridge Quantum's tket.
Each pass must preserve the statevector modulo global phase:

    RemoveRedundancies, CliffordSimp, FullPeepholeOptimise, SynthesiseTket

A random circuit is built, one pass applied, and the state compared
before/after with the N1 amplitude oracle. Divergences shrink to a
1-minimal witness and go to the bug zoo.

Usage:
    python3 tools/fuzz_tket_passes.py --iterations 4000 --seed 42
"""

import argparse
import json
import pathlib
import random
import sys

import numpy as np
from pytket import Circuit
from pytket.passes import (
    CliffordSimp,
    FullPeepholeOptimise,
    RemoveRedundancies,
    SynthesiseTket,
)

TOL = 1e-9

PASSES = {
    "RemoveRedundancies": RemoveRedundancies,
    "CliffordSimp": CliffordSimp,
    "FullPeepholeOptimise": FullPeepholeOptimise,
    "SynthesiseTket": SynthesiseTket,
}

# Matched to fuzz_pytket.py (half-turn params for tket rotations).
GATES = [
    ("Rz", 1, 1, 5), ("U1", 1, 1, 5), ("CRz", 2, 1, 5), ("CU1", 2, 1, 5),
    ("S", 1, 0, 5), ("Sdg", 1, 0, 5), ("T", 1, 0, 5), ("Tdg", 1, 0, 5),
    ("SX", 1, 0, 5), ("SXdg", 1, 0, 5),
    ("H", 1, 0, 4),
    ("Rx", 1, 1, 2), ("Ry", 1, 1, 2),
    ("CX", 2, 0, 3), ("CY", 2, 0, 1), ("CZ", 2, 0, 1),
    ("CSX", 2, 0, 2), ("CRx", 2, 1, 1), ("CRy", 2, 1, 1),
    ("X", 1, 0, 1), ("Y", 1, 0, 1), ("Z", 1, 0, 1),
    ("SWAP", 2, 0, 1), ("CCX", 3, 0, 1),
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


def divergence(ops, num_qubits, pass_factory):
    try:
        c = build(num_qubits, ops)
        sv = c.get_statevector()
        c2 = c.copy()
        pass_factory().apply(c2)
        sv2 = c2.get_statevector()
        if sv.shape != sv2.shape:
            return float("inf")
        return n1_amplitude(sv, sv2)
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
    for name, qubits, params in ops:
        p = f"({','.join(f'{x:.4f}' for x in params)})" if params else ""
        out.append(f"{name}{p} {qubits}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=4000)
    ap.add_argument("--max-qubits", type=int, default=4)
    ap.add_argument("--max-depth", type=int, default=12)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--pass", dest="pass_name", choices=list(PASSES) + ["all"], default="all")
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    args = ap.parse_args()

    import pytket

    which = list(PASSES) if args.pass_name == "all" else [args.pass_name]
    print(f"fuzz_tket_passes — pytket {pytket.__version__}")
    print(f"  passes: {which}")
    print(f"  iterations: {args.iterations}  seed: {args.seed}  tol: {args.tol:g}\n")

    total = 0
    for pname in which:
        factory = PASSES[pname]
        found = 0
        worst = 0.0
        for i in range(args.iterations):
            rng = random.Random(args.seed + i)
            nq = rng.randint(1, args.max_qubits)
            depth = rng.randint(1, args.max_depth)
            ops = random_ops(rng, nq, depth)

            d = divergence(ops, nq, factory)
            if d is None:
                continue
            worst = max(worst, d)
            if d <= args.tol:
                continue

            found += 1
            total += 1
            minimal = shrink(ops, nq, lambda c: (divergence(c, nq, factory) or 0.0) > args.tol)
            final = divergence(minimal, nq, factory) or 0.0
            print(f"── [{pname}] divergence #{found} circuit {i}, N1 = {d:.3e}")
            print(f"  minimal witness ({len(minimal)} gates):")
            for line in describe(minimal):
                print(f"    {line}")

            entry = {
                "id": f"tket-{pname}-{found:03}",
                "target": f"tket {pname} pass",
                "pytket_version": pytket.__version__,
                "generator_seed": args.seed + i,
                "num_qubits": nq,
                "gates": [{"gate": g, "qubits": q, "params": p} for g, q, p in minimal],
                "n1_amplitude_distance": final,
            }
            args.zoo_dir.mkdir(exist_ok=True)
            out = args.zoo_dir / f"{entry['id']}.json"
            out.write_text(json.dumps(entry, indent=2))
            print(f"  zoo entry: {out}\n")

        tag = "CLEAN" if worst <= args.tol else "BUG"
        print(f"[{pname}] {found} divergence(s), worst N1 = {worst:.3e}  ({tag})\n")

    print(f"Total divergences across passes: {total}")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
