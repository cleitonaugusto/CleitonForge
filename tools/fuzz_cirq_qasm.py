#!/usr/bin/env python3
"""Differential fuzzer for Cirq's OpenQASM 2 export.

The property is the same one that found the CUDA-Q bug: the OpenQASM 2 that
`cirq.qasm` emits must implement the same operator as the circuit Cirq itself
simulates.

Why this surface and not Cirq's optimizers: those were already fuzzed from this
repo, 3000 circuits, and came back clean. That result covered half the code.
Export was never touched, and export is where the last two bugs of this kind
lived, in CUDA-Q and in Lift. Translators get read and assumed correct because
"it only prints"; optimizers get fuzzed because everyone expects bugs there.

Cirq's gate set is wider than OpenQASM 2's, so most of these have to be
decomposed or rewritten on the way out. That is the part worth testing.

Two halves, because Cirq lives in its own venv here and the oracle stack lives
in the system Python. Neither half trusts the other's toolchain:

    ~/.venvs/cirq/bin/python tools/fuzz_cirq_qasm.py --emit  cases.jsonl
    python3 tools/fuzz_cirq_qasm.py --check cases.jsonl
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import random

# Exponents a compiler special-cases. Cirq expresses most gates as a power, so
# the exponent is where its rounding and special-casing decisions live: 1.0 is
# the plain gate, 0.5 is a square root, and 2.0 is the identity up to phase.
SPECIAL_EXP = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, -0.5, -1.0, 4.0]
SPECIAL_ANGLE = sorted({s * g for s in (1, -1) for g in
                        (0.0, math.pi / 8, math.pi / 4, math.pi / 2,
                         math.pi, 3 * math.pi / 2, 2 * math.pi, 4 * math.pi)})

# (name, num_qubits, kind) where kind says what the parameter means.
GATES = [
    ("H", 1, None), ("X", 1, None), ("Y", 1, None), ("Z", 1, None),
    ("S", 1, None), ("T", 1, None),
    ("XPow", 1, "exp"), ("YPow", 1, "exp"), ("ZPow", 1, "exp"),
    ("Rx", 1, "angle"), ("Ry", 1, "angle"), ("Rz", 1, "angle"),
    ("PhasedXPow", 1, "two_exp"),
    ("CNOT", 2, None), ("CZ", 2, None), ("SWAP", 2, None), ("ISWAP", 2, None),
    ("CZPow", 2, "exp"), ("SwapPow", 2, "exp"), ("ISwapPow", 2, "exp"),
    ("CCX", 3, None), ("CCZ", 3, None), ("CSWAP", 3, None),
]


def draw(rng: random.Random, kind: str, boundary_prob: float):
    if kind == "angle":
        if boundary_prob and rng.random() < boundary_prob:
            return [rng.choice(SPECIAL_ANGLE)]
        return [rng.uniform(-math.pi, math.pi)]
    if boundary_prob and rng.random() < boundary_prob:
        pick = lambda: rng.choice(SPECIAL_EXP)  # noqa: E731
    else:
        pick = lambda: rng.uniform(-2.0, 2.0)  # noqa: E731
    return [pick(), pick()] if kind == "two_exp" else [pick()]


def random_ops(rng: random.Random, num_qubits: int, depth: int,
               boundary_prob: float) -> list:
    eligible = [g for g in GATES if g[1] <= num_qubits]
    ops = []
    for _ in range(depth):
        name, nq, kind = rng.choice(eligible)
        qubits = rng.sample(range(num_qubits), nq)
        params = draw(rng, kind, boundary_prob) if kind else []
        ops.append((name, qubits, params))
    return ops


def build_circuit(ops, num_qubits: int):
    import cirq

    q = cirq.LineQubit.range(num_qubits)
    table = {
        "H": lambda p: cirq.H, "X": lambda p: cirq.X, "Y": lambda p: cirq.Y,
        "Z": lambda p: cirq.Z, "S": lambda p: cirq.S, "T": lambda p: cirq.T,
        "XPow": lambda p: cirq.XPowGate(exponent=p[0]),
        "YPow": lambda p: cirq.YPowGate(exponent=p[0]),
        "ZPow": lambda p: cirq.ZPowGate(exponent=p[0]),
        "Rx": lambda p: cirq.rx(p[0]), "Ry": lambda p: cirq.ry(p[0]),
        "Rz": lambda p: cirq.rz(p[0]),
        "PhasedXPow": lambda p: cirq.PhasedXPowGate(phase_exponent=p[0], exponent=p[1]),
        "CNOT": lambda p: cirq.CNOT, "CZ": lambda p: cirq.CZ,
        "SWAP": lambda p: cirq.SWAP, "ISWAP": lambda p: cirq.ISWAP,
        "CZPow": lambda p: cirq.CZPowGate(exponent=p[0]),
        "SwapPow": lambda p: cirq.SwapPowGate(exponent=p[0]),
        "ISwapPow": lambda p: cirq.ISwapPowGate(exponent=p[0]),
        "CCX": lambda p: cirq.CCX, "CCZ": lambda p: cirq.CCZ,
        "CSWAP": lambda p: cirq.CSWAP,
    }
    circuit = cirq.Circuit()
    for name, qubits, params in ops:
        circuit.append(table[name](params).on(*[q[i] for i in qubits]))
    return circuit, q


def emit(args) -> int:
    import cirq
    import numpy as np

    def record_for(ops, num_qubits: int) -> dict:
        rec = {"ops": ops, "num_qubits": num_qubits}
        try:
            circuit, order = build_circuit(ops, num_qubits)
            sv = cirq.Simulator().simulate(circuit, qubit_order=order).final_state_vector
            sv = np.asarray(sv)
            if len(sv) != 2 ** num_qubits:
                raise ValueError(f"state has {len(sv)}, expected {2 ** num_qubits}")
            rec["state"] = [[float(a.real), float(a.imag)] for a in sv]
            rec["qasm2"] = cirq.qasm(circuit)
        except Exception as e:  # noqa: BLE001 — an export failure is a finding
            rec["error"] = f"{type(e).__name__}: {e}"
        return rec

    out = pathlib.Path(args.emit)
    written = errors = 0
    with out.open("w") as fh:
        # Calibration first, and it is Cirq's own answer rather than a literal:
        # X on qubit 0 of two is asymmetric under bit reversal, so it is the
        # case that pins the convention. Recording it here means the check half
        # compares two toolchains instead of comparing Qiskit to itself.
        cal = record_for([("X", [0], [])], 2)
        cal["calibration"] = True
        fh.write(json.dumps(cal) + "\n")

        for i in range(args.iterations):
            rng = random.Random(args.seed + i)
            depth = rng.randint(args.min_depth, args.max_depth)
            ops = random_ops(rng, args.qubits, depth, args.boundary_prob)
            rec = record_for(ops, args.qubits)
            rec.update({"i": i, "seed": args.seed + i, "depth": depth,
                        "min_depth": args.min_depth, "max_depth": args.max_depth,
                        "boundary_prob": args.boundary_prob})
            if "error" in rec:
                errors += 1
            else:
                written += 1
            fh.write(json.dumps(rec) + "\n")
    print(f"emitted {written} cases, {errors} export failures -> {out}")
    return 0


EXACT_THRESHOLD = 1e-12
MISMATCH_THRESHOLD = 1e-6


def qubit_map(qasm: str) -> list | None:
    """Which original qubit each QASM register slot stands for.

    Cirq drops qubits an circuit never touches and renumbers what is left, so a
    3 qubit circuit that only uses q(0) and q(2) exports as `qreg q[2]` with
    q(2) written as q[1]. It records the mapping in a comment:

        // Qubits: [q(0), q(2)]

    That is a documented convention, not a defect, and comparing the state
    without honouring it invents divergences on any circuit with an idle qubit.
    """
    import re

    m = re.search(r"^// Qubits: \[(.*)\]", qasm, re.M)
    if not m:
        return None
    return [int(x) for x in re.findall(r"q\((\d+)\)", m.group(1))]


def embed(sub, mapping: list, num_qubits: int):
    """Lift a state over the exported qubits into the full register.

    Qubits Cirq dropped were never acted on, so they stay in |0>. Bit b of the
    exported index belongs to original qubit mapping[b]; both sides count from
    the most significant bit, which is Cirq's convention.
    """
    import numpy as np

    k = len(mapping)
    full = np.zeros(2 ** num_qubits, dtype=complex)
    for j, amp in enumerate(sub):
        if amp == 0:
            continue
        idx = 0
        for b in range(k):
            if (j >> (k - 1 - b)) & 1:
                idx |= 1 << (num_qubits - 1 - mapping[b])
        full[idx] = amp
    return full


def rebuild(qasm: str, num_qubits: int):
    """Operator the emitted QASM implements, in Cirq's bit order and width.

    Parse and simulate are kept apart so an emitter convention (a trailing
    measurement, say) is not reported as hundreds of parse bugs.
    """
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector

    qc = QuantumCircuit.from_qasm_str(qasm)
    # Cirq puts qubit 0 in the most significant bit, Qiskit in the least.
    # Measured, not assumed: X on qubit 0 of two lands on index 2 under Cirq.
    sub = Statevector.from_instruction(qc).reverse_qargs().data
    mapping = qubit_map(qasm)
    if mapping is None or len(mapping) != qc.num_qubits:
        raise ValueError("could not read the // Qubits: mapping from the export")
    return qc, embed(sub, mapping, num_qubits)


def divergence(want, got) -> float:
    """Largest amplitude difference modulo global phase.

    Not fidelity: `1 - |<a|b>|` is quadratic in the error, so a 1e-6 bound on it
    only fires above an angle of about 1.4e-3.
    """
    import numpy as np

    overlap = np.vdot(want, got)
    if abs(overlap) < 1e-12:
        return float(np.max(np.abs(want - got)))
    phase = overlap / abs(overlap)
    return float(np.max(np.abs(want - got * np.conj(phase))))


def check(args) -> int:
    import numpy as np

    records = [json.loads(l) for l in pathlib.Path(args.check).read_text().splitlines()]

    cal = next((r for r in records if r.get("calibration")), None)
    if cal is None or "state" not in cal:
        print("NO CALIBRATION CASE. Re-run --emit. Refusing to report verdicts.")
        return 2
    try:
        _, cal_got = rebuild(cal["qasm2"], cal["num_qubits"])
    except Exception as e:  # noqa: BLE001
        print(f"CALIBRATION UNPARSABLE: {type(e).__name__}: {e}")
        return 2
    cal_want = np.array([complex(r, im) for r, im in cal["state"]])
    if divergence(cal_want, cal_got) > 1e-9:
        print("CALIBRATION FAILED: Cirq and the rebuilt circuit disagree on X "
              "applied to qubit 0, the case that distinguishes the two bit "
              "orders. Every verdict below would be meaningless.")
        return 2

    mismatches = crashes = export_failures = checked = imprecise = 0
    for case in records:
        if case.get("calibration"):
            continue
        if "error" in case:
            export_failures += 1
            print(f"  [export failed] case {case['i']}: {case['error'][:140]}")
            continue
        try:
            qc, got = rebuild(case["qasm2"], case["num_qubits"])
        except Exception as e:  # noqa: BLE001
            crashes += 1
            print(f"  [rebuild failed] case {case['i']}: {type(e).__name__}: {e}")
            continue
        want = np.array([complex(r, im) for r, im in case["state"]])
        if len(got) != len(want):
            crashes += 1
            print(f"  [width mismatch] case {case['i']}: {len(want)} vs {len(got)}")
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

    print(f"\nchecked {checked}, {mismatches} mismatch, {imprecise} lossy, "
          f"{crashes} unusable, {export_failures} export failures")
    if not checked:
        print("NOTHING WAS CHECKED. A run that tested nothing is not a pass.")
    return 1 if (mismatches or crashes or export_failures or not checked) else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit", metavar="FILE")
    ap.add_argument("--check", metavar="FILE")
    ap.add_argument("--iterations", type=int, default=300)
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
    raise SystemExit(main())
