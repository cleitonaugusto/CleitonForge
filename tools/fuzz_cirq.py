#!/usr/bin/env python3
"""Metamorphic fuzzing of Cirq's optimizing transformers — Phase 3.

Same architecture as the other harnesses, targeting Google's Cirq. Each
transformer must preserve the statevector modulo global phase:

    eject_z, eject_phased_paulis, merge_single_qubit_gates_to_phxz,
    optimize_for_target_gateset (CZ target)

The phase-ejecting transformers push Z/phase gates through the circuit,
which is exactly the kind of pass where a rotation-sign convention slip
would surface. A random circuit is built, one transformer applied, and
the state compared before/after with the N1 amplitude oracle.

Note: a stub for mpl_toolkits.mplot3d is installed before importing cirq
to work around a system/pip matplotlib split on this machine (cirq only
needs it for Bloch-sphere plotting, which we don't use).

Usage:
    python3 tools/fuzz_cirq.py --iterations 4000 --seed 42
"""

import argparse
import json
import pathlib
import random
import sys
import types

import numpy as np

# Work around a broken system mpl_toolkits shadowing pip matplotlib 3.10.
_m = types.ModuleType("mpl_toolkits.mplot3d")
_ax = types.ModuleType("mpl_toolkits.mplot3d.axes3d")
_ax.Axes3D = object
_m.axes3d = _ax
sys.modules["mpl_toolkits.mplot3d"] = _m
sys.modules["mpl_toolkits.mplot3d.axes3d"] = _ax

import cirq  # noqa: E402

TOL = 1e-9

TRANSFORMERS = {
    "eject_z": lambda c: cirq.eject_z(c),
    "eject_phased_paulis": lambda c: cirq.eject_phased_paulis(c),
    "merge_single_qubit_gates_to_phxz": lambda c: cirq.merge_single_qubit_gates_to_phxz(c),
    "optimize_for_cz_target": lambda c: cirq.optimize_for_target_gateset(
        c, gateset=cirq.CZTargetGateset()
    ),
}


def gate_pool():
    # (name, builder(params)->gate, n_qubits, n_params, weight). sx/sxdg via
    # XPowGate(±0.5) to mirror the Qiskit bug's gate family.
    return [
        ("rz", lambda p: cirq.rz(p[0]), 1, 1, 5),
        ("phase", lambda p: cirq.ZPowGate(exponent=p[0] / np.pi), 1, 1, 5),
        ("s", lambda p: cirq.S, 1, 0, 5),
        ("s_dag", lambda p: cirq.S**-1, 1, 0, 5),
        ("t", lambda p: cirq.T, 1, 0, 5),
        ("t_dag", lambda p: cirq.T**-1, 1, 0, 5),
        ("sx", lambda p: cirq.XPowGate(exponent=0.5), 1, 0, 5),
        ("sx_dag", lambda p: cirq.XPowGate(exponent=-0.5), 1, 0, 5),
        ("h", lambda p: cirq.H, 1, 0, 4),
        ("rx", lambda p: cirq.rx(p[0]), 1, 1, 2),
        ("ry", lambda p: cirq.ry(p[0]), 1, 1, 2),
        ("x", lambda p: cirq.X, 1, 0, 1),
        ("y", lambda p: cirq.Y, 1, 0, 1),
        ("z", lambda p: cirq.Z, 1, 0, 1),
        ("cnot", lambda p: cirq.CNOT, 2, 0, 3),
        ("cz", lambda p: cirq.CZ, 2, 0, 1),
        ("swap", lambda p: cirq.SWAP, 2, 0, 1),
        ("ccx", lambda p: cirq.CCX, 3, 0, 1),
    ]


POOL = gate_pool()


def random_ops(rng, num_qubits, depth):
    eligible = [g for g in POOL if g[2] <= num_qubits]
    weights = [g[4] for g in eligible]
    ops = []
    for _ in range(depth):
        name, builder, nq, npar, _ = rng.choices(eligible, weights=weights)[0]
        qubits = rng.sample(range(num_qubits), nq)
        params = [rng.uniform(-np.pi, np.pi) for _ in range(npar)]
        ops.append((name, builder, qubits, params))
    return ops


def build(num_qubits, ops):
    q = cirq.LineQubit.range(num_qubits)
    c = cirq.Circuit()
    for _name, builder, qubits, params in ops:
        c.append(builder(params).on(*[q[i] for i in qubits]))
    return c, q


def state(num_qubits, ops):
    c, q = build(num_qubits, ops)
    return cirq.Simulator(dtype=np.complex128).simulate(c, qubit_order=q).final_state_vector


def state_of_circuit(c, q):
    return cirq.Simulator(dtype=np.complex128).simulate(c, qubit_order=q).final_state_vector


def n1_amplitude(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if a.shape != b.shape:
        return float("inf")
    overlap = np.vdot(a, b)
    phase = overlap / abs(overlap) if abs(overlap) > 1e-15 else 1.0
    return float(np.max(np.abs(a * phase - b)))


def divergence(num_qubits, ops, transformer):
    try:
        c, q = build(num_qubits, ops)
        sv0 = state_of_circuit(c, q)
        c2 = transformer(c)
        sv1 = state_of_circuit(c2, q)
        return n1_amplitude(sv0, sv1)
    except Exception as e:  # noqa: BLE001
        print(f"  [crash] {type(e).__name__}: {e}", file=sys.stderr)
        return None


def shrink(num_qubits, ops, still_fails):
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
    for name, _b, qubits, params in ops:
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
    ap.add_argument("--transformer", choices=list(TRANSFORMERS) + ["all"], default="all")
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    args = ap.parse_args()

    which = list(TRANSFORMERS) if args.transformer == "all" else [args.transformer]
    print(f"fuzz_cirq — cirq {cirq.__version__}")
    print(f"  transformers: {which}")
    print(f"  iterations: {args.iterations}  seed: {args.seed}  tol: {args.tol:g}\n")

    total = 0
    for tname in which:
        transformer = TRANSFORMERS[tname]
        found = 0
        worst = 0.0
        for i in range(args.iterations):
            rng = random.Random(args.seed + i)
            nq = rng.randint(1, args.max_qubits)
            depth = rng.randint(1, args.max_depth)
            ops = random_ops(rng, nq, depth)

            d = divergence(nq, ops, transformer)
            if d is None:
                continue
            worst = max(worst, d)
            if d <= args.tol:
                continue

            found += 1
            total += 1
            minimal = shrink(nq, ops, lambda c: (divergence(nq, c, transformer) or 0.0) > args.tol)
            final = divergence(nq, minimal, transformer) or 0.0
            print(f"── [{tname}] divergence #{found} circuit {i}, N1 = {d:.3e}")
            print(f"  minimal witness ({len(minimal)} gates):")
            for line in describe(minimal):
                print(f"    {line}")

            entry = {
                "id": f"cirq-{tname}-{found:03}",
                "target": f"Cirq {tname} transformer",
                "cirq_version": cirq.__version__,
                "generator_seed": args.seed + i,
                "num_qubits": nq,
                "gates": [
                    {"gate": n, "qubits": q, "params": p} for n, _b, q, p in minimal
                ],
                "n1_amplitude_distance": final,
            }
            args.zoo_dir.mkdir(exist_ok=True)
            out = args.zoo_dir / f"{entry['id']}.json"
            out.write_text(json.dumps(entry, indent=2))
            print(f"  zoo entry: {out}\n")

        tag = "CLEAN" if worst <= args.tol else "BUG"
        print(f"[{tname}] {found} divergence(s), worst N1 = {worst:.3e}  ({tag})\n")

    print(f"Total divergences across transformers: {total}")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
