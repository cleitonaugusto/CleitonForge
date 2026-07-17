#!/usr/bin/env python3
"""QCEC-backed equivalence oracle.

Replaces the hand-rolled amplitude oracle (max |a*phase - b| over the full
2^n x 2^n unitary). That oracle capped campaigns at 4 qubits: at 12 it takes
28s per circuit, and beyond that it is infeasible. This one is flat to 20
qubits (~0.03s), which is what makes layout and routing passes reachable.

Configuration is not a matter of taste. Measured over 120 random transpiled
circuits:

    default                                 120/120 agreement, NON-deterministic
    run_simulation_checker=False, seed=N    120/120 agreement, deterministic, 2.4x faster

The simulation checker is probabilistic. It returns "probably equivalent",
and when it disagrees with the ZX checker the manager reports no_information.
Turning it off removes both the nondeterminism and that whole verdict class.

Phase semantics (verified, do not trust the type stub). The .pyi docstrings for
equivalent_up_to_phase and equivalent_up_to_global_phase are swapped: the latter
is documented as "global or relative", which would make this oracle blind to
exactly the bug class we hunt. The behaviour is correct. Checked directly:

    H vs H*e^{i1.234}          -> equivalent                    accepted
    Rz(pi/2) vs P(pi/2)        -> equivalent_up_to_global_phase accepted  (a real global phase)
    Rz(pi/2) vs Rz(-pi/2)      -> not_equivalent                flagged
    S vs Sdg                   -> not_equivalent                flagged
    I vs Z                     -> not_equivalent                flagged
    Rz(1.0) vs Rz(-1.0)        -> not_equivalent                flagged  (the quantrs2 bug)

Global phase is quotiented, relative phase is not. That is the property this
whole project depends on.
"""

from __future__ import annotations

import enum

from mqt import qcec

DEFAULT_SEED = 42
DEFAULT_TIMEOUT = 60.0


class Verdict(enum.Enum):
    """Three outcomes. UNKNOWN is not a bug -- keep them apart."""

    EQUIVALENT = "equivalent"
    BUG = "bug"
    UNKNOWN = "unknown"


# EquivalenceCriterion -> Verdict. Names taken from the enum, not from the
# docstrings, which are wrong (see module docstring).
_ACCEPT = {"equivalent", "equivalent_up_to_phase", "equivalent_up_to_global_phase"}
_FLAG = {"not_equivalent"}
# no_information: timeout or nothing ran.
# probably_equivalent: simulation checker gave up (off by default here).
# probably_not_equivalent: the ZX checker could not reduce to identity. That is
# a failure to prove, NOT a proof of failure -- counting it as a bug would
# manufacture false positives.
_UNKNOWN = {"no_information", "probably_equivalent", "probably_not_equivalent"}


def criterion_name(result) -> str:
    return str(result.equivalence).replace("EquivalenceCriterion.", "")


def verify(qc, other, seed: int = DEFAULT_SEED, timeout: float = DEFAULT_TIMEOUT):
    """Run the check and return the raw QCEC result."""
    return qcec.verify(
        qc,
        other,
        run_simulation_checker=False,
        seed=seed,
        timeout=timeout,
    )


def verdict(qc, other, seed: int = DEFAULT_SEED, timeout: float = DEFAULT_TIMEOUT) -> Verdict:
    """EQUIVALENT / BUG / UNKNOWN for a pair of circuits.

    Never raises for an inconclusive check: that is UNKNOWN. It does propagate
    exceptions from QCEC itself, since a crash in the checker is its own finding.
    """
    name = criterion_name(verify(qc, other, seed=seed, timeout=timeout))
    if name in _ACCEPT:
        return Verdict.EQUIVALENT
    if name in _FLAG:
        return Verdict.BUG
    if name in _UNKNOWN:
        return Verdict.UNKNOWN
    raise RuntimeError(f"unhandled EquivalenceCriterion: {name!r}")


def is_bug(qc, other, seed: int = DEFAULT_SEED, timeout: float = DEFAULT_TIMEOUT) -> bool:
    """Convenience for the shrinker, which only cares whether it still fails."""
    return verdict(qc, other, seed=seed, timeout=timeout) is Verdict.BUG


# Widths where building the full 2^n x 2^n unitary is affordable. This is the
# oracle we replaced, kept only as a second opinion on the small circuits a
# shrinker produces, never as the primary check.
CONFIRM_MAX_QUBITS = 10


def confirm_bug(qc, other) -> bool | None:
    """Second opinion on a BUG verdict. True if real, False if QCEC is wrong.

    Returns None when the circuit is too wide to check this way.

    QCEC is the better oracle but it is not infallible, and it fails in a way
    that matters here: comparing a parametrised gate against the same operator
    packed into a UnitaryGate, it reports not_equivalent for small angles even
    though the two are bit-identical. Measured on ConsolidateBlocks output,
    crz(theta) against its own matrix:

        theta <= 1e-4   ->  not_equivalent, and the true error is exactly 0.0
        theta >= 1e-3   ->  equivalent

    Its decision diagram appears to prune entries below a tolerance, which
    collapses a nearly-identity matrix while the rotation gate keeps its angle.
    That produced a finding that looked new and was not there at all.

    A witness is small by the time it is worth reporting, so the exact
    comparison always fits. Cross-checking costs nothing and this project has no
    business trusting a single oracle -- that is the whole argument it makes to
    everyone else.
    """
    if qc.num_qubits > CONFIRM_MAX_QUBITS:
        return None
    from qiskit.quantum_info import Operator

    return not Operator(qc).equiv(Operator(other))


if __name__ == "__main__":
    import math

    from qiskit import QuantumCircuit

    cases = []

    a = QuantumCircuit(1); a.h(0)
    b = QuantumCircuit(1); b.h(0); b.global_phase = 1.234
    cases.append(("H vs H*e^{i1.234}", a, b, Verdict.EQUIVALENT))

    a = QuantumCircuit(1); a.rz(math.pi / 2, 0)
    b = QuantumCircuit(1); b.p(math.pi / 2, 0)
    cases.append(("Rz(pi/2) vs P(pi/2)", a, b, Verdict.EQUIVALENT))

    a = QuantumCircuit(1); a.rz(1.0, 0)
    b = QuantumCircuit(1); b.rz(-1.0, 0)
    cases.append(("Rz(1.0) vs Rz(-1.0)", a, b, Verdict.BUG))

    a = QuantumCircuit(1); a.s(0)
    b = QuantumCircuit(1); b.sdg(0)
    cases.append(("S vs Sdg", a, b, Verdict.BUG))

    a = QuantumCircuit(2); a.sxdg(1); a.sx(0); a.sxdg(1)
    from qiskit import transpile
    b = transpile(a, basis_gates=["rz", "sx", "x", "cx"], optimization_level=3,
                  approximation_degree=1.0, seed_transpiler=7)
    cases.append(("qiskit #16594 witness", a, b, Verdict.BUG))

    bad = 0
    for label, x, y, want in cases:
        got = verdict(x, y)
        ok = got is want
        bad += not ok
        print(f"  {'ok  ' if ok else 'FAIL'} {label:<26} {got.value:<10} (want {want.value})")
    raise SystemExit(1 if bad else 0)
