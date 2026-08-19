#!/usr/bin/env python3
"""Differential fuzzer for CUDA-Q's compilation to OpenQASM 2.

The property: the OpenQASM 2 that `cudaq.translate` emits must implement the
same operator as the kernel CUDA-Q itself simulates. Anything else means the
compiler and the simulator disagree about what the program is, and only one of
them can be right.

That path is worth testing because it is not a printer: reaching OpenQASM 2 runs
the decomposition patterns in lib/Optimizer/Transforms/DecompositionPatterns.cpp.
Decomposition tables are where the last two bugs of this kind lived: a sign on
the first rotation of an RX sequence, and a measurement conjugation whose
closing half was never removed.

What this covers, stated plainly because it is smaller than it sounds: the gate
pool below is single-control at most, has no adjoint variants, and no exp_pauli.
Nothing here counts which decomposition patterns were actually reached, so any
claim about pattern coverage would be unmeasured.

Runs in two halves because CUDA-Q needs Python 3.11 and our oracle stack lives
on 3.10:

    ~/.venvs/cudaq/bin/python tools/fuzz_cudaq.py --emit  cases.jsonl
    python3 tools/fuzz_cudaq.py --check cases.jsonl

The emit half generates kernels, records the state CUDA-Q computes, and asks
for the OpenQASM 2. The check half rebuilds the operator from that QASM and
compares. Neither half trusts the other's toolchain.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import random
import sys

# Gates reachable from the CUDA-Q Python DSL, as (name, num_qubits, num_params).
#
# An earlier version carried a fourth column claiming to mark which gates have
# no OpenQASM 2 equivalent. It was wrong (qelib1.inc does provide swap, cy, cz,
# ch, ccx, crz and cp) and nothing read it, so it is gone rather than fixed:
# a field that lies and is never consulted is worse than no field.
GATES = [
    ("h", 1, 0),
    ("x", 1, 0),
    ("y", 1, 0),
    ("z", 1, 0),
    ("s", 1, 0),
    ("t", 1, 0),
    ("rx", 1, 1),
    ("ry", 1, 1),
    ("rz", 1, 1),
    ("r1", 1, 1),
    ("swap", 2, 0),
    ("x.ctrl", 2, 0),
    ("y.ctrl", 2, 0),
    ("z.ctrl", 2, 0),
    ("h.ctrl", 2, 0),
    ("s.ctrl", 2, 0),
    ("t.ctrl", 2, 0),
    ("rx.ctrl", 2, 1),
    ("ry.ctrl", 2, 1),
    ("rz.ctrl", 2, 1),
    ("r1.ctrl", 2, 1),
    ("swap.ctrl", 3, 0),
]

# Angles a compiler special-cases, same rationale as tools/generator.py: a
# uniform draw hits an exact multiple of pi/2 with probability zero, and that
# is where rotation bugs live.
SPECIAL = [0.0, math.pi / 8, math.pi / 4, math.pi / 2, 3 * math.pi / 4,
           math.pi, 3 * math.pi / 2, 2 * math.pi, 4 * math.pi]
SPECIAL = sorted({s * sign for s in SPECIAL for sign in (1, -1)})


def random_angle(rng: random.Random, boundary_prob: float) -> float:
    if boundary_prob and rng.random() < boundary_prob:
        return rng.choice(SPECIAL)
    return rng.uniform(-math.pi, math.pi)


def random_ops(rng: random.Random, num_qubits: int, depth: int,
               boundary_prob: float) -> list:
    eligible = [g for g in GATES if g[1] <= num_qubits]
    ops = []
    for _ in range(depth):
        name, nq, npar = rng.choice(eligible)
        qubits = rng.sample(range(num_qubits), nq)
        params = [random_angle(rng, boundary_prob) for _ in range(npar)]
        ops.append((name, qubits, params))
    return ops


def kernel_source(name: str, ops, num_qubits: int) -> str:
    """CUDA-Q kernel source. Angles are baked in as literals, since translating
    a kernel with arguments to OpenQASM 2 is not supported."""
    body = [f"    q = cudaq.qvector({num_qubits})"]
    for gate, qubits, params in ops:
        args = [repr(float(p)) for p in params] + [f"q[{i}]" for i in qubits]
        body.append(f"    {gate}({', '.join(args)})")
    return f"@cudaq.kernel\ndef {name}():\n" + "\n".join(body) + "\n"


def emit(args) -> int:
    """Generate kernels, record what CUDA-Q simulates, ask for the OpenQASM 2.

    The kernels are written to a real module rather than exec'd: the
    @cudaq.kernel decorator reads the function's source text back off disk, so
    a string compiled in memory is rejected.
    """
    import importlib.util
    import tempfile

    import cudaq

    tmpdir = pathlib.Path(tempfile.mkdtemp(prefix="cudaq_fuzz_"))

    def build_kernel(name: str, ops, num_qubits: int):
        """Import one kernel in its own module.

        @cudaq.kernel bridges the AST to MLIR at decoration time, so a kernel
        that trips a front-end bug raises during import. Importing them all from
        one module means the first such kernel kills the run and takes every
        other result with it, which loses exactly the case worth keeping.
        """
        src = "import cudaq\nfrom cudaq import *\n\n" + kernel_source(name, ops, num_qubits)
        path = tmpdir / f"{name}.py"
        path.write_text(src)
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return getattr(mod, name)

    def record_for(name: str, ops, num_qubits: int) -> dict:
        rec = {"ops": ops, "num_qubits": num_qubits}
        try:
            kernel = build_kernel(name, ops, num_qubits)
            state = cudaq.get_state(kernel)
            width = 2 ** num_qubits
            amps = [complex(state[j]) for j in range(width)]
            if len(state) != width:
                raise ValueError(f"state has {len(state)} amplitudes, expected {width}")
            rec["state"] = [[a.real, a.imag] for a in amps]
            rec["qasm2"] = cudaq.translate(kernel, format="openqasm2")
        except Exception as e:  # noqa: BLE001 — a compiler crash is a finding
            rec["error"] = f"{type(e).__name__}: {e}"
        return rec

    out = pathlib.Path(args.emit)
    written = errors = 0
    with out.open("w") as fh:
        # Calibration case first, and it is CUDA-Q's own answer rather than a
        # literal: x on qubit 0 of two, which is asymmetric under bit reversal.
        # The check half compares against these recorded amplitudes, so the
        # guard fails if either side ever changes convention.
        cal = record_for("calib", [("x", [0], [])], 2)
        cal["calibration"] = True
        fh.write(json.dumps(cal) + "\n")

        for i in range(args.iterations):
            rng = random.Random(args.seed + i)
            depth = rng.randint(args.min_depth, args.max_depth)
            ops = random_ops(rng, args.qubits, depth, args.boundary_prob)
            rec = record_for(f"gen_{i}", ops, args.qubits)
            rec.update({"i": i, "seed": args.seed + i, "depth": depth,
                        "min_depth": args.min_depth, "max_depth": args.max_depth,
                        "boundary_prob": args.boundary_prob})
            if "error" in rec:
                errors += 1
            else:
                written += 1
            fh.write(json.dumps(rec) + "\n")
    print(f"emitted {written} cases, {errors} compiler errors -> {out}")
    print(f"kernels kept under {tmpdir}")
    return 0


def rebuild(qasm: str):
    """Operator that the emitted OpenQASM 2 implements, in CUDA-Q's bit order.

    Parse and simulate are separated on purpose: if the emitter ever appends a
    measurement or a reset, the QASM parses fine and only the simulation fails,
    and calling both under one `except` would report 200 emitter-convention
    failures as 200 parse bugs.
    """
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector

    qc = QuantumCircuit.from_qasm_str(qasm)
    # CUDA-Q puts qubit 0 in the most significant bit, Qiskit in the least.
    # Measured, not assumed: x(q[0]) on two qubits lands on index 2 under
    # CUDA-Q and index 1 under Qiskit.
    return qc, Statevector.from_instruction(qc).reverse_qargs().data


# Two thresholds, because there are two things to separate and reporting them
# as one number hides both.
#
# The emitter prints rotation angles with about seven significant digits: a
# kernel asking for rx(-1.2397473823778045) comes back as rx(-1.239747e+00).
# That alone moves amplitudes by up to ~1e-7, so anything below MISMATCH is
# precision loss rather than a wrong circuit, and is counted separately instead
# of being either hidden or called a bug.
#
# Measured on a 60-circuit run: divergences cluster at machine epsilon (<5e-16)
# or at 2.6e-9 and above, with nothing in between.
EXACT_THRESHOLD = 1e-12
MISMATCH_THRESHOLD = 1e-6


def divergence(want, got) -> float:
    """Largest amplitude difference modulo global phase.

    Not fidelity. `1 - |<a|b>|` is quadratic in the error, so a 1e-6 threshold
    on it only fires above an angle of about 1.4e-3, which would score a
    decomposition that is off by a milliradian as a pass. This is the same trap
    tools/oracle_qcec.py documents for Operator.equiv's default rtol.
    """
    import numpy as np

    overlap = np.vdot(want, got)
    if abs(overlap) < 1e-12:
        return float(np.max(np.abs(want - got)))
    # Conjugate, do not multiply. If got = e^{i.phi} want then vdot(want, got)
    # is e^{i.phi}, so aligning means got * conj(phase). Multiplying instead
    # applies the rotation twice and scores two identical states as maximally
    # different: [1,0] against [i,0] came out as 2.0, the largest value the
    # metric can take.
    phase = overlap / abs(overlap)
    return float(np.max(np.abs(want - got * np.conj(phase))))


def check(args) -> int:
    import numpy as np

    path = pathlib.Path(args.check)
    records = [json.loads(l) for l in path.read_text().splitlines()]

    # The calibration case is CUDA-Q's own answer for an asymmetric circuit, so
    # this compares the two toolchains rather than comparing Qiskit to itself.
    # An earlier version hardcoded the expected amplitudes here, which made the
    # guard tautological: it could not fail, and deleting the bit reversal in
    # the comparison below left it passing.
    cal = next((r for r in records if r.get("calibration")), None)
    if cal is None or "state" not in cal:
        print("NO CALIBRATION CASE in the input. Re-run --emit with a current "
              "version. Refusing to report verdicts that cannot be trusted.")
        return 2
    cal_want = np.array([complex(r, im) for r, im in cal["state"]])
    try:
        _, cal_got = rebuild(cal["qasm2"])
    except Exception as e:  # noqa: BLE001
        print(f"CALIBRATION UNPARSABLE: {type(e).__name__}: {e}")
        return 2
    if divergence(cal_want, cal_got) > 1e-9:
        print("CALIBRATION FAILED: CUDA-Q and the rebuilt circuit disagree on "
              "x applied to qubit 0, which is the case that distinguishes the "
              "two bit orders. Every verdict below would be meaningless.")
        return 2
    mismatches = crashes = compiler_errors = checked = imprecise = 0
    for case in records:
        if case.get("calibration"):
            continue
        if "error" in case:
            compiler_errors += 1
            print(f"  [compiler error] case {case['i']}: {case['error']}")
            continue
        try:
            qc, got = rebuild(case["qasm2"])
        except Exception as e:  # noqa: BLE001
            crashes += 1
            print(f"  [rebuild failed] case {case['i']}: {type(e).__name__}: {e}")
            continue

        want = np.array([complex(r, im) for r, im in case["state"]])
        if len(got) != len(want):
            crashes += 1
            print(f"  [width mismatch] case {case['i']}: CUDA-Q reported "
                  f"{len(want)} amplitudes, the emitted QASM has {len(got)} "
                  f"({qc.num_qubits} qubits). Not scored as a divergence.")
            continue
        checked += 1

        d = divergence(want, got)
        if d > MISMATCH_THRESHOLD:
            mismatches += 1
            print(f"\n  [MISMATCH] case {case['i']} (seed {case['seed']}), "
                  f"max amplitude difference {d:.3e}")
            for name, qubits, params in case["ops"]:
                ps = ", ".join(f"{p:.6f}" for p in params)
                print(f"      {name}({ps}) on {qubits}")
        elif d > EXACT_THRESHOLD:
            imprecise += 1

    print(f"\nchecked {checked}, {mismatches} mismatch, "
          f"{imprecise} within angle-precision loss, "
          f"{crashes} unusable, {compiler_errors} compiler errors")
    if not checked:
        print("NOTHING WAS CHECKED. A run that tested nothing is not a pass.")
    # compiler_errors counts too: the emit half calls a crash a finding, and a
    # green exit on a run where every kernel failed to compile is the worst
    # failure mode a fuzzer has.
    return 1 if (mismatches or crashes or compiler_errors or not checked) else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit", metavar="FILE", help="run under the cudaq venv")
    ap.add_argument("--check", metavar="FILE", help="run under the oracle venv")
    ap.add_argument("--iterations", type=int, default=200)
    ap.add_argument("--qubits", type=int, default=3)
    ap.add_argument("--min-depth", type=int, default=2)
    ap.add_argument("--max-depth", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--boundary-prob", type=float, default=0.5)
    args = ap.parse_args()
    if args.emit:
        return emit(args)
    if args.check:
        return check(args)
    ap.error("pass --emit or --check")
    return 2


if __name__ == "__main__":
    sys.exit(main())
