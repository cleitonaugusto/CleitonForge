#!/usr/bin/env python3
"""Wide metamorphic fuzzing of the Qiskit transpiler, QCEC-backed.

Every campaign so far ran at 4 qubits, all-to-all, because the amplitude
oracle had to build the full 2^n x 2^n unitary (28s at 12 qubits, infeasible
above that). At that width the interesting passes never fire: with no coupling
map there is nothing to route, and layout is trivial. We found #16594 through
that keyhole.

This harness targets the passes we have never tested. It transpiles against a
real heavy-hex coupling map, so SabreLayout, routing and swap insertion all
engage, and it checks the result with QCEC, which is flat to 20 qubits.

The property is unchanged: transpilation must preserve the unitary up to global
phase. Qiskit tracks global_phase, and routing is a permutation QCEC accounts
for, so any not_equivalent verdict is a correctness bug.

approximation_degree=1.0 stays. It is a declared transformation: without it the
transpiler is licensed to drop sub-threshold gates and QCEC would rightly call
that inequivalent, flooding the run with findings that are contract, not bug.

Usage:
    python3 tools/fuzz_qiskit_wide.py --iterations 500 --lattice 3 --seed 42
"""

from __future__ import annotations

import argparse
import json
import pathlib
import random
import sys
import time

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from generator import build, random_ops, shrink, to_qasm2  # noqa: E402
from oracle_qcec import Verdict, criterion_name, verdict, verify  # noqa: E402

IBM_BASIS = ["rz", "sx", "x", "cx"]


def heavy_hex(distance: int) -> CouplingMap:
    """Full heavy-hex lattice. distance 3 -> 19 qubits, 5 -> 57.

    The circuit width is the lattice size, not a free parameter. Truncating to
    a prefix does not work: the induced subgraph of heavy-hex on 0..15 is
    disconnected, so a 16-qubit circuit cannot be routed on it. Taking the
    whole lattice also keeps the original and transpiled circuits the same
    width, which is what lets QCEC compare them without ancilla handling.
    """
    cm = CouplingMap.from_heavy_hex(distance)
    if not cm.is_connected():
        raise ValueError(f"heavy-hex distance {distance} is not connected")
    return cm


def compile_pair(ops, num_qubits: int, coupling: CouplingMap, opt_level: int):
    """(original, transpiled) or None if the circuit cannot be compiled."""
    qc = build(num_qubits, ops)
    tqc = transpile(
        qc,
        basis_gates=IBM_BASIS,
        coupling_map=coupling,
        optimization_level=opt_level,
        approximation_degree=1.0,
        seed_transpiler=7,
    )
    return qc, tqc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=500)
    ap.add_argument("--lattice", type=int, default=3,
                    help="heavy-hex distance: 3 -> 19 qubits, 5 -> 57")
    ap.add_argument("--qubits", type=int, default=None,
                    help="width for --all-to-all; ignored with a coupling map "
                         "(the lattice fixes it)")
    ap.add_argument("--min-depth", type=int, default=20)
    ap.add_argument("--max-depth", type=int, default=60)
    ap.add_argument("--timeout", type=float, default=10.0,
                    help="per-check seconds. Measured at 19 qubits: depth<=60 is "
                         "0.18s/circuit and always conclusive; at depth 100 the "
                         "decision diagrams blow up (31s/circuit, 2/8 timing out)")
    ap.add_argument("--opt-level", type=int, default=3)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    ap.add_argument("--all-to-all", action="store_true",
                    help="no coupling map (the old narrow regime)")
    args = ap.parse_args()

    import qiskit

    if args.all_to_all:
        coupling = None
        num_qubits = args.qubits or 16
        topo = "all-to-all"
    else:
        coupling = heavy_hex(args.lattice)
        num_qubits = coupling.size()
        topo = f"heavy-hex d={args.lattice} ({len(coupling.get_edges())} edges)"

    print(f"fuzz_qiskit_wide — qiskit {qiskit.__version__}, oracle: QCEC")
    print(f"  {args.iterations} circuits, {num_qubits} qubits, depth "
          f"{args.min_depth}-{args.max_depth}, O{args.opt_level}, {topo}")
    print(f"  seed {args.seed}\n")

    bugs = unknowns = errors = 0
    t_start = time.perf_counter()

    for i in range(args.iterations):
        # Before the verdict checks: most iterations end in a `continue`, so a
        # progress report placed after them never runs on a clean campaign.
        if i and i % 50 == 0:
            el = time.perf_counter() - t_start
            print(f"[{i}/{args.iterations}] {el:.0f}s, {el / i:.2f}s/circuit, "
                  f"{bugs} bug(s), {unknowns} unknown, {errors} error(s)")

        rng = random.Random(args.seed + i)
        depth = rng.randint(args.min_depth, args.max_depth)
        ops = random_ops(rng, num_qubits, depth)

        try:
            qc, tqc = compile_pair(ops, num_qubits, coupling, args.opt_level)
        except Exception as e:  # noqa: BLE001 — a transpiler crash is a finding
            errors += 1
            print(f"  [crash] circuit {i}: {type(e).__name__}: {e}")
            continue

        try:
            v = verdict(qc, tqc, timeout=args.timeout)
        except Exception as e:  # noqa: BLE001
            errors += 1
            print(f"  [oracle crash] circuit {i}: {type(e).__name__}: {e}")
            continue

        if v is Verdict.EQUIVALENT:
            continue
        if v is Verdict.UNKNOWN:
            unknowns += 1
            print(f"  [unknown] circuit {i} — inconclusive, not counted as a bug")
            continue

        bugs += 1
        print(f"\n── divergence #{bugs} at circuit {i} ({depth} gates)")

        def still_fails(cand):
            try:
                pair = compile_pair(cand, num_qubits, coupling, args.opt_level)
            except Exception:  # noqa: BLE001
                return False
            try:
                return verdict(*pair, timeout=args.timeout) is Verdict.BUG
            except Exception:  # noqa: BLE001
                return False

        minimal = shrink(ops, still_fails)
        print(f"  minimal witness ({len(minimal)} gates):")
        for line in to_qasm2(minimal, num_qubits).splitlines()[3:]:
            print(f"    {line}")

        qc_m, tqc_m = compile_pair(minimal, num_qubits, coupling, args.opt_level)
        res = verify(qc_m, tqc_m)
        entry = {
            "id": f"qiskit-wide-{num_qubits}q-{bugs:03}",
            "qiskit_version": qiskit.__version__,
            "oracle": "qcec (run_simulation_checker=False)",
            "criterion": criterion_name(res),
            "topology": topo,
            "optimization_level": args.opt_level,
            "generator_seed": args.seed + i,
            "num_qubits": num_qubits,
            "gates": [{"gate": g, "qubits": qs, "params": ps} for g, qs, ps in minimal],
            "qasm2": to_qasm2(minimal, num_qubits),
        }
        args.zoo_dir.mkdir(exist_ok=True)
        out = args.zoo_dir / f"{entry['id']}.json"
        out.write_text(json.dumps(entry, indent=2))
        print(f"  zoo entry: {out}\n")

    el = time.perf_counter() - t_start
    print(f"\nResult: {bugs} divergence(s), {unknowns} inconclusive, "
          f"{errors} error(s) in {args.iterations} circuits ({el:.0f}s).")
    return 1 if bugs else 0


if __name__ == "__main__":
    sys.exit(main())
