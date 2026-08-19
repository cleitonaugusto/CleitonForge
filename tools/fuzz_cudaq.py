#!/usr/bin/env python3
"""Differential fuzzer for CUDA-Q's compilation to OpenQASM 2.

The property: the OpenQASM 2 that `cudaq.translate` emits must implement the
same operator as the kernel CUDA-Q itself simulates. Anything else means the
compiler and the simulator disagree about what the program is, and only one of
them can be right.

That path is worth testing because it is not a printer. OpenQASM 2 has no r1,
no exp_pauli, no multi-controlled anything, so getting there runs the
decomposition patterns in lib/Optimizer/Transforms/DecompositionPatterns.cpp,
66 of them. Decomposition tables are where the last two bugs of this kind lived:
a sign on the first rotation of an RX sequence, and a measurement conjugation
whose closing half was never removed.

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

# Gates reachable from the CUDA-Q Python DSL. `qasm2` is the name emitted when
# the gate survives to OpenQASM 2 unchanged; None means it has to be
# decomposed, which is the interesting case.
#   (name, num_qubits, num_params, needs_decomposition)
GATES = [
    ("h", 1, 0, False),
    ("x", 1, 0, False),
    ("y", 1, 0, False),
    ("z", 1, 0, False),
    ("s", 1, 0, False),
    ("t", 1, 0, False),
    ("rx", 1, 1, False),
    ("ry", 1, 1, False),
    ("rz", 1, 1, False),
    ("r1", 1, 1, True),       # phase gate, no QASM2 equivalent
    ("swap", 2, 0, True),
    ("x.ctrl", 2, 0, False),  # cx
    ("y.ctrl", 2, 0, True),
    ("z.ctrl", 2, 0, True),
    ("h.ctrl", 2, 0, True),
    ("s.ctrl", 2, 0, True),
    ("t.ctrl", 2, 0, True),
    ("rx.ctrl", 2, 1, True),
    ("ry.ctrl", 2, 1, True),
    ("rz.ctrl", 2, 1, True),
    ("r1.ctrl", 2, 1, True),
    ("swap.ctrl", 3, 0, True),
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
        name, nq, npar, _ = rng.choice(eligible)
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

    cases = []
    for i in range(args.iterations):
        rng = random.Random(args.seed + i)
        depth = rng.randint(args.min_depth, args.max_depth)
        cases.append((i, random_ops(rng, args.qubits, depth, args.boundary_prob)))

    header = ["import cudaq", "from cudaq import *", ""]
    src = "\n".join(header + [kernel_source(f"gen_{i}", ops, args.qubits)
                              for i, ops in cases])

    tmpdir = pathlib.Path(tempfile.mkdtemp(prefix="cudaq_fuzz_"))
    modpath = tmpdir / "generated_kernels.py"
    modpath.write_text(src)
    spec = importlib.util.spec_from_file_location("generated_kernels", modpath)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["generated_kernels"] = mod
    spec.loader.exec_module(mod)

    out = pathlib.Path(args.emit)
    written = errors = 0
    with out.open("w") as fh:
        for i, ops in cases:
            record = {"i": i, "seed": args.seed + i, "ops": ops,
                      "num_qubits": args.qubits}
            try:
                kernel = getattr(mod, f"gen_{i}")
                state = cudaq.get_state(kernel)
                record["state"] = [[complex(state[j]).real, complex(state[j]).imag]
                                   for j in range(2 ** args.qubits)]
                record["qasm2"] = cudaq.translate(kernel, format="openqasm2")
                written += 1
            except Exception as e:  # noqa: BLE001 — a compiler crash is a finding
                errors += 1
                record["error"] = f"{type(e).__name__}: {e}"
            fh.write(json.dumps(record) + "\n")
    print(f"emitted {written} cases, {errors} compiler errors -> {out}")
    print(f"kernels kept at {modpath}")
    return 0


def calibrate() -> bool:
    """Prove the comparison can tell qubit 0 from qubit 1 before trusting it.

    The first version of this harness compared the two state vectors directly
    and reported 3 mismatches in 5 circuits. The bug was mine: CUDA-Q indexes
    qubit 0 as the most significant bit and Qiskit as the least. The smoke test
    that missed it used a Bell pair, which is invariant under bit reversal, so
    it passed under either convention and proved nothing.

    This runs an asymmetric circuit whose answer differs between the two
    conventions, so a run cannot start unless the reversal is right.
    """
    import numpy as np
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector

    qc = QuantumCircuit(2)
    qc.x(0)
    correct = Statevector.from_instruction(qc).reverse_qargs().data
    naive = Statevector.from_instruction(qc).data
    # CUDA-Q's answer for x on qubit 0 of two: index 2, measured.
    cudaq_answer = np.array([0, 0, 1, 0], dtype=complex)
    return (abs(abs(np.vdot(cudaq_answer, correct)) - 1.0) < 1e-9
            and abs(np.vdot(cudaq_answer, naive)) < 1e-9)


def check(args) -> int:
    import numpy as np
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector

    if not calibrate():
        print("CALIBRATION FAILED: the comparison cannot distinguish qubit "
              "order, so every verdict below would be meaningless. Aborting.")
        return 2

    path = pathlib.Path(args.check)
    mismatches = crashes = compiler_errors = checked = 0
    for line in path.read_text().splitlines():
        case = json.loads(line)
        if "error" in case:
            compiler_errors += 1
            print(f"  [compiler error] case {case['i']}: {case['error']}")
            continue
        try:
            qc = QuantumCircuit.from_qasm_str(case["qasm2"])
            # CUDA-Q puts qubit 0 in the most significant bit, Qiskit in the
            # least. Measured, not assumed: x(q[0]) on two qubits lands on
            # index 2 under CUDA-Q and index 1 under Qiskit. Comparing without
            # this reversal reports a mismatch on any circuit that is not
            # symmetric under bit reversal, which a Bell pair happens to be.
            got = Statevector.from_instruction(qc).reverse_qargs().data
        except Exception as e:  # noqa: BLE001
            crashes += 1
            print(f"  [unparsable qasm] case {case['i']}: {type(e).__name__}: {e}")
            continue

        want = np.array([complex(r, im) for r, im in case["state"]])
        checked += 1

        # Compare up to global phase, which no measurement can see.
        overlap = np.vdot(want, got)
        if abs(abs(overlap) - 1.0) > 1e-6:
            mismatches += 1
            print(f"\n  [MISMATCH] case {case['i']} (seed {case['seed']}), "
                  f"|<want|got>| = {abs(overlap):.6f}")
            for name, qubits, params in case["ops"]:
                ps = ", ".join(f"{p:.6f}" for p in params)
                print(f"      {name}({ps}) on {qubits}")

    print(f"\nchecked {checked}, {mismatches} mismatch, "
          f"{crashes} unparsable, {compiler_errors} compiler errors")
    return 1 if (mismatches or crashes) else 0


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
