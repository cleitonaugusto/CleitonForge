#!/usr/bin/env python3
"""Phase-sensitive differential fuzzing of Qiskit — Phase 3 campaign.

Two experiments, mirroring the cforge-fuzz oracle hierarchy:

1. TRANSPILER (metamorphic): |psi(C)> vs |psi(transpile(C, O3))>.
   The transpiler must preserve the unitary exactly (Qiskit tracks
   global_phase), so any N1 divergence modulo global phase is a
   correctness bug in a >5k-star project.

2. AER (differential): exact Statevector(C) vs AerSimulator
   statevector method on the same circuit.

Divergences are shrunk to 1-minimal witnesses (greedy gate removal)
and written to the bug zoo as JSON.

Usage:
    python3 tools/fuzz_qiskit.py --iterations 2000 --seed 42
"""

import argparse
import json
import math
import pathlib
import random
import sys

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator

# Same weighted pool as cforge-fuzz generator.rs: phase-family and
# discrete complex gates 5x, so bugged-rotation x correct-complex-gate
# mixtures dominate the search.
GATES = [
    # (name, n_qubits, n_params, weight)
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

TOL = 1e-9  # transpiler must be exact up to f64 accumulation


def random_ops(rng: random.Random, num_qubits: int, depth: int):
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


def build(num_qubits: int, ops) -> QuantumCircuit:
    qc = QuantumCircuit(num_qubits)
    for name, qubits, params in ops:
        getattr(qc, name)(*params, *qubits)
    return qc


# ── Oracles (same math as cforge-fuzz/src/oracle.rs) ────────────────

def n1_amplitude(a: np.ndarray, b: np.ndarray) -> float:
    """Max amplitude distance modulo one global phase."""
    overlap = np.vdot(a, b)
    phase = overlap / abs(overlap) if abs(overlap) > 1e-15 else 1.0
    return float(np.max(np.abs(a * phase - b)))


def n2_probability(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(np.abs(a) ** 2 - np.abs(b) ** 2)))


# ── Experiments ──────────────────────────────────────────────────────

IBM_BASIS = ["rz", "sx", "x", "cx"]


def transpiler_divergence(ops, num_qubits: int) -> float | None:
    """N1 distance between original and O3-transpiled circuit, or None
    if the circuit fails to build/transpile.

    The transpiler is *allowed* to permute qubits (ElidePermutations at
    O2+ removes swaps and records the permutation in the layout), so the
    oracle must undo the declared permutation before comparing —
    Operator.from_circuit does exactly that. Comparing raw statevectors
    instead yields ~4% false positives (measured: 119/3000 circuits,
    every witness containing a swap).
    """
    try:
        from qiskit.quantum_info import Operator

        qc = build(num_qubits, ops)
        op_orig = Operator(qc).data
        tqc = transpile(
            qc, basis_gates=IBM_BASIS, optimization_level=3, seed_transpiler=7
        )
        # from_circuit applies the transpiled layout/routing permutation.
        op_trans = Operator.from_circuit(tqc).data
        # N1 on full unitaries, modulo one global phase.
        return n1_amplitude(op_orig.reshape(-1), op_trans.reshape(-1))
    except Exception as e:  # noqa: BLE001 — any crash is also a finding
        print(f"  [crash] {type(e).__name__}: {e}", file=sys.stderr)
        return None


def aer_divergence(ops, num_qubits: int) -> float | None:
    """N1 distance between exact Statevector and Aer statevector."""
    try:
        qc = build(num_qubits, ops)
        sv_exact = np.asarray(Statevector(qc))
        qc2 = qc.copy()
        qc2.save_statevector()
        backend = AerSimulator(method="statevector")
        # Light translation to Aer's native gate set (no optimization),
        # since a few stdgates (e.g. ch) have no direct Aer instruction.
        qc2 = transpile(qc2, backend, optimization_level=0)
        result = backend.run(qc2, shots=1).result()
        sv_aer = np.asarray(result.get_statevector())
        return n1_amplitude(sv_exact, sv_aer)
    except Exception as e:  # noqa: BLE001
        print(f"  [crash] {type(e).__name__}: {e}", file=sys.stderr)
        return None


def shrink(ops, num_qubits: int, still_fails) -> list:
    """Greedy 1-minimal reduction, same as cforge-fuzz shrinker pass 2."""
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


def to_qasm2(ops, num_qubits: int) -> str:
    lines = ['OPENQASM 2.0;', 'include "qelib1.inc";', f"qreg q[{num_qubits}];"]
    for name, qubits, params in ops:
        p = f"({','.join(f'{x:.12f}' for x in params)})" if params else ""
        q = ",".join(f"q[{i}]" for i in qubits)
        lines.append(f"{name}{p} {q};")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=2000)
    ap.add_argument("--max-qubits", type=int, default=4)
    ap.add_argument("--max-depth", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    ap.add_argument(
        "--experiment", choices=["transpiler", "aer", "both"], default="both"
    )
    args = ap.parse_args()

    import qiskit
    import qiskit_aer

    experiments = []
    if args.experiment in ("transpiler", "both"):
        experiments.append(("qiskit-transpiler-O3", transpiler_divergence))
    if args.experiment in ("aer", "both"):
        experiments.append(("qiskit-aer-statevector", aer_divergence))

    print(f"fuzz_qiskit — qiskit {qiskit.__version__}, aer {qiskit_aer.__version__}")
    print(f"  iterations: {args.iterations}  seed: {args.seed}  tol: {args.tol:g}")
    print(f"  experiments: {[name for name, _ in experiments]}\n")

    found = 0
    worst: dict[str, float] = {name: 0.0 for name, _ in experiments}

    for i in range(args.iterations):
        rng = random.Random(args.seed + i)
        nq = rng.randint(1, args.max_qubits)
        depth = rng.randint(1, args.max_depth)
        ops = random_ops(rng, nq, depth)

        for name, experiment in experiments:
            dist = experiment(ops, nq)
            if dist is None:
                continue
            worst[name] = max(worst[name], dist)
            if dist <= args.tol:
                continue

            found += 1
            minimal = shrink(
                ops, nq,
                lambda cand: (experiment(cand, nq) or 0.0) > args.tol,
            )
            final = experiment(minimal, nq) or 0.0
            print(f"── divergence #{found} [{name}] circuit {i}, N1 = {dist:.3e}")
            print(f"  minimal witness ({len(minimal)} gates):")
            for line in to_qasm2(minimal, nq).splitlines()[3:]:
                print(f"    {line}")

            entry = {
                "id": f"{name}-{found:03}",
                "qiskit_version": qiskit.__version__,
                "aer_version": qiskit_aer.__version__,
                "experiment": name,
                "generator_seed": args.seed + i,
                "num_qubits": nq,
                "gates": [
                    {"gate": g, "qubits": qs, "params": ps}
                    for g, qs, ps in minimal
                ],
                "qasm2": to_qasm2(minimal, nq),
                "n1_amplitude_distance": final,
            }
            args.zoo_dir.mkdir(exist_ok=True)
            out = args.zoo_dir / f"{entry['id']}.json"
            out.write_text(json.dumps(entry, indent=2))
            print(f"  zoo entry: {out}\n")

        if (i + 1) % 500 == 0:
            print(f"[{i + 1}/{args.iterations}] worst N1 so far: "
                  + ", ".join(f"{k}={v:.2e}" for k, v in worst.items()))

    print(f"\nResult: {found} divergence(s) in {args.iterations} circuits.")
    print("Worst N1 distance per experiment (below tol = clean):")
    for k, v in worst.items():
        print(f"  {k}: {v:.3e}")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
