#!/usr/bin/env python3
"""Fuzz individual transpiler passes, QCEC-backed.

The wide campaign transpiles end to end, which is a blunt instrument: a
divergence tells you the pipeline is wrong, and the shrinker then has to work
out which of ~30 passes did it. #16594 lived in one pass, and testing that pass
alone would have pointed straight at it.

So this runs one pass at a time. The property is the one every optimization
pass owes: applying it must preserve the circuit's unitary up to global phase.
Qiskit tracks global_phase, so a not_equivalent verdict is a correctness bug in
that pass, and the pass name is the diagnosis -- no shrinking needed to locate
it, only to minimize the witness.

Passes here take no required arguments and are pure rewrites, so nothing is
licensed to change except gate structure. That is what makes the oracle sound
without the quotients the end-to-end harness needs (no layout, no routing, no
approximate synthesis). The exception is RemoveIdentityEquivalent, which is
explicitly approximate: it is pinned to approximation_degree=1.0 for the same
reason transpile() is in the wide harness.

CommutativeCancellation is included as a positive control, not as a target: it
carries a known bug (#16594), so a run that reports it clean means the harness
is broken -- but only once the run is long enough for that to mean anything.
The fault surfaces in about 0.7% of random 4-qubit circuits, so 60 iterations
expect 0.4 hits and a clean control there is the likely outcome, not evidence
of anything. The check waits until the expected count is at least 3.

Usage:
    python3 tools/fuzz_qiskit_passes.py --iterations 2000
    python3 tools/fuzz_qiskit_passes.py --pass InverseCancellation --iterations 5000
"""

from __future__ import annotations

import argparse
import json
import pathlib
import random
import sys
import time

from qiskit.transpiler import PassManager
from qiskit.transpiler import passes as qpasses

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from generator import build, random_ops, shrink, to_qasm2  # noqa: E402
from oracle_qcec import Verdict, criterion_name, verdict, verify  # noqa: E402


def _pass_factories():
    """name -> zero-arg constructor. Only pure rewrites belong here."""
    f = {}
    if hasattr(qpasses, "InverseCancellation"):
        f["InverseCancellation"] = lambda: qpasses.InverseCancellation()
    if hasattr(qpasses, "CommutativeCancellation"):
        # Positive control: known bug #16594.
        f["CommutativeCancellation"] = lambda: qpasses.CommutativeCancellation()
    if hasattr(qpasses, "CommutativeInverseCancellation"):
        f["CommutativeInverseCancellation"] = lambda: qpasses.CommutativeInverseCancellation()
        f["CommutativeInverseCancellation-matrix"] = (
            lambda: qpasses.CommutativeInverseCancellation(matrix_based=True))
    if hasattr(qpasses, "Optimize1qGatesDecomposition"):
        f["Optimize1qGatesDecomposition"] = lambda: qpasses.Optimize1qGatesDecomposition(
            basis=["rz", "sx", "x"])
    if hasattr(qpasses, "Optimize1qGatesSimpleCommutation"):
        f["Optimize1qGatesSimpleCommutation"] = lambda: qpasses.Optimize1qGatesSimpleCommutation(
            basis=["rz", "sx", "x"])
    if hasattr(qpasses, "ConsolidateBlocks"):
        f["ConsolidateBlocks"] = lambda: qpasses.ConsolidateBlocks(force_consolidate=True)
    if hasattr(qpasses, "OptimizeCliffordT"):
        f["OptimizeCliffordT"] = lambda: qpasses.OptimizeCliffordT()
    if hasattr(qpasses, "RemoveIdentityEquivalent"):
        # Approximate by contract; pin it off, as the wide harness does.
        f["RemoveIdentityEquivalent"] = lambda: qpasses.RemoveIdentityEquivalent(
            approximation_degree=1.0)
    return f


PASSES = _pass_factories()


def run_pass(factory, ops, num_qubits):
    qc = build(num_qubits, ops)
    out = PassManager([factory()]).run(qc)
    return qc, out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass", dest="pass_name", choices=sorted(PASSES), default=None,
                    help="default: every pass")
    ap.add_argument("--iterations", type=int, default=2000)
    ap.add_argument("--qubits", type=int, default=4)
    ap.add_argument("--min-depth", type=int, default=4)
    ap.add_argument("--max-depth", type=int, default=14)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--timeout", type=float, default=10.0)
    ap.add_argument("--boundary-prob", type=float, default=0.0,
                    help="probability of drawing an angle from the special "
                         "values a compiler branches on, rather than uniform. "
                         "A uniform draw hits an exact multiple of pi/2 with "
                         "probability zero (measured: 0 in 200000), which is "
                         "where rotation bugs live")
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--zoo-dir", type=pathlib.Path, default=pathlib.Path("bug-zoo"))
    args = ap.parse_args()

    import qiskit

    targets = [args.pass_name] if args.pass_name else sorted(PASSES)
    excluded = frozenset(args.exclude)

    print(f"fuzz_qiskit_passes — qiskit {qiskit.__version__}, oracle: QCEC")
    print(f"  {args.iterations} circuits/pass, {args.qubits} qubits, "
          f"depth {args.min_depth}-{args.max_depth}, seed {args.seed}")
    if excluded:
        print(f"  excluding: {sorted(excluded)}")
    print(f"  passes: {targets}\n")

    totals = {}
    for name in targets:
        factory = PASSES[name]
        bugs = unknowns = errors = 0
        first_witness = None
        t0 = time.perf_counter()

        for i in range(args.iterations):
            rng = random.Random(args.seed + i)
            depth = rng.randint(args.min_depth, args.max_depth)
            ops = random_ops(rng, args.qubits, depth, exclude=excluded, boundary_prob=args.boundary_prob)
            try:
                qc, out = run_pass(factory, ops, args.qubits)
            except Exception as e:  # noqa: BLE001 — a pass crash is a finding
                errors += 1
                if errors <= 3:
                    print(f"  [{name}] crash on circuit {i}: {type(e).__name__}: {e}")
                continue
            try:
                v = verdict(qc, out, timeout=args.timeout)
            except Exception as e:  # noqa: BLE001
                errors += 1
                if errors <= 3:
                    print(f"  [{name}] oracle crash on circuit {i}: {type(e).__name__}: {e}")
                continue

            if v is Verdict.EQUIVALENT:
                continue
            if v is Verdict.UNKNOWN:
                unknowns += 1
                continue

            bugs += 1
            if first_witness is None:
                def still_fails(cand):
                    try:
                        return verdict(*run_pass(factory, cand, args.qubits),
                                       timeout=args.timeout) is Verdict.BUG
                    except Exception:  # noqa: BLE001
                        return False

                minimal = shrink(ops, still_fails)
                first_witness = minimal
                print(f"\n── [{name}] divergence at circuit {i}, "
                      f"minimal witness ({len(minimal)} gates):")
                for line in to_qasm2(minimal, args.qubits).splitlines()[3:]:
                    print(f"    {line}")
                qc_m, out_m = run_pass(factory, minimal, args.qubits)
                entry = {
                    "id": f"pass-{name}-{args.seed}",
                    "qiskit_version": qiskit.__version__,
                    "oracle": "qcec (run_simulation_checker=False)",
                    "pass": name,
                    "criterion": criterion_name(verify(qc_m, out_m)),
                    "generator_seed": args.seed + i,
                    "num_qubits": args.qubits,
                    "gates": [{"gate": g, "qubits": q, "params": p} for g, q, p in minimal],
                    "qasm2": to_qasm2(minimal, args.qubits),
                }
                args.zoo_dir.mkdir(exist_ok=True)
                out_path = args.zoo_dir / f"{entry['id']}.json"
                out_path.write_text(json.dumps(entry, indent=2))
                print(f"    zoo: {out_path}\n")

        el = time.perf_counter() - t0
        totals[name] = (bugs, unknowns, errors)
        rate = 100.0 * bugs / max(args.iterations, 1)
        print(f"{name:<38} {bugs:>4} bug ({rate:5.2f}%)  "
              f"{unknowns:>3} unknown  {errors:>3} err  {el:>6.0f}s")

    print("\n== summary ==")
    control = totals.get("CommutativeCancellation")
    for name, (b, u, e) in totals.items():
        tag = ""
        if name == "CommutativeCancellation":
            tag = "  (positive control: #16594)"
        print(f"  {name:<38} {b:>4} bug  {u:>3} unknown  {e:>3} err{tag}")

    # Measured reachability of #16594 at 4 qubits, from the earlier campaigns.
    CONTROL_RATE = 0.007
    expected = args.iterations * CONTROL_RATE
    if control is not None and control[0] == 0:
        if expected >= 3:
            print(f"\nWARNING: the positive control found nothing, and at "
                  f"{args.iterations} iterations it should have found about "
                  f"{expected:.0f}. #16594 is known to live in "
                  f"CommutativeCancellation on this version, so this harness is "
                  f"not testing what it thinks it is. Do not trust the other rows.")
            return 2
        print(f"\nNote: the positive control found nothing, but at "
              f"{args.iterations} iterations it only expects {expected:.1f}, so "
              f"that says nothing either way. Run at least "
              f"{int(3 / CONTROL_RATE)} to make the control meaningful.")
    return 1 if any(b for b, _, _ in totals.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
