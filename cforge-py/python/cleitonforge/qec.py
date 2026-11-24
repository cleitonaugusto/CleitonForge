"""Check whether a QEC compiler pass changed the circuit it was given.

A pass that rewrites an error-correction circuit is usually accepted when the
detector error model (DEM) before and after are equal. That test cannot see a
whole class of wrong rewrites: a deterministic Pauli inserted by the pass moves
the measurement outcomes and the detector reference frame together, so the DEM
comes out identical while the circuit computes something else. The smallest
case is four instructions long, see `witness()`.

This module checks the circuit itself, through its stabiliser flows, and reports
the DEM comparison next to it so the blind spot is visible instead of silent:

    from cleitonforge import qec

    verdict = qec.check_pass(before, after)
    verdict.equivalent   # the flows of the two circuits match, both ways
    verdict.dem_equal    # what a DEM comparison would have said
    verdict.dem_blind    # DEM equal although the circuit changed
    verdict.witness      # a flow one circuit has and the other lacks

It also measures how much a given check can see, by applying mutations that are
known to change a circuit and counting which ones the check reports:

    qec.oracle_power(qec.dem_differs, [circuit])

Circuits are `stim.Circuit`, or any object with an `as_stim_circuit()` method,
which covers `deltakit_circuit.Circuit`.

Requires `stim` (`pip install cleitonforge[qec]`).
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Callable, Iterable, Optional, Tuple

import stim

__all__ = [
    "MUTATIONS",
    "Verdict",
    "check_pass",
    "dem_differs",
    "dem_equal",
    "flows_differ",
    "flows_equivalent",
    "oracle_power",
    "to_stim",
    "witness",
]


def to_stim(circuit) -> stim.Circuit:
    """Return `circuit` as a `stim.Circuit`, converting it when it knows how."""
    if isinstance(circuit, stim.Circuit):
        return circuit
    convert = getattr(circuit, "as_stim_circuit", None)
    if convert is None:
        raise TypeError(
            f"expected a stim.Circuit or an object with as_stim_circuit(), "
            f"got {type(circuit).__name__}"
        )
    return convert()


def _closed(circuit: stim.Circuit, num_qubits: int, discard_final_state: bool) -> stim.Circuit:
    closed = circuit.without_noise().flattened()
    if discard_final_state:
        # The state a circuit leaves behind is usually thrown away. Resetting
        # every qubit at the end says so, and keeps a pass that removes a gate
        # nobody reads from being reported as a change.
        closed.append("R", range(num_qubits))
    return closed


def flows_equivalent(
    before, after, *, discard_final_state: bool = True
) -> Tuple[bool, Optional[stim.Flow]]:
    """Compare two circuits by their stabiliser flows, in both directions.

    Noise is ignored: this checks what the circuit computes, not how it fails.

    Args:
        before: The circuit given to the pass.
        after: The circuit the pass returned.
        discard_final_state: Treat the state left on the qubits at the end as
            discarded. True suits a memory experiment, where only measurement
            results matter. False also compares the output state.

    Returns:
        `(True, None)` when every flow of each circuit is a flow of the other.
        Otherwise `(False, flow)`, with a flow that one of them has and the
        other lacks.
    """
    a, b = to_stim(before), to_stim(after)
    n = max(a.num_qubits, b.num_qubits)
    a, b = _closed(a, n, discard_final_state), _closed(b, n, discard_final_state)
    if a.num_measurements != b.num_measurements:
        return False, None
    for x, y in ((a, b), (b, a)):
        for flow in x.flow_generators():
            if not y.has_flow(flow):
                return False, flow
    return True, None


def flows_differ(before, after) -> bool:
    """Oracle form of `flows_equivalent`: True when the circuits differ."""
    return not flows_equivalent(before, after)[0]


def _mechanisms(dem: stim.DetectorErrorModel) -> dict:
    """Error mechanisms of a DEM, keyed by the detectors and observables they flip.

    Detector coordinates and the order of the lines are left out, since a
    decoder does not depend on them. Mechanisms with the same symptoms are
    merged into one, with the probability that an odd number of them fires.
    """
    merged = {}
    for ins in dem.flattened():
        if ins.type != "error":
            continue
        key = frozenset(str(t) for t in ins.targets_copy())
        p = ins.args_copy()[0]
        q = merged.get(key, 0.0)
        merged[key] = q + p - 2 * q * p
    return merged


def dem_equal(before, after, *, tolerance: float = 1e-9) -> bool:
    """Whether the two circuits have the same detector error model, as a decoder
    sees it: the same detectors, observables and error mechanisms.

    This is the comparison commonly used to accept a pass. It is reported here
    because `check_pass` needs to show when it would have been fooled.
    A circuit whose DEM cannot be built counts as different.
    """
    try:
        dem_a = to_stim(before).detector_error_model(decompose_errors=False)
        dem_b = to_stim(after).detector_error_model(decompose_errors=False)
    except ValueError:
        return False
    if (dem_a.num_detectors, dem_a.num_observables) != (dem_b.num_detectors, dem_b.num_observables):
        return False
    a, b = _mechanisms(dem_a), _mechanisms(dem_b)
    return a.keys() == b.keys() and all(abs(a[k] - b[k]) <= tolerance for k in a)


def dem_differs(before, after) -> bool:
    """Oracle form of `dem_equal`: True when the DEMs differ."""
    return not dem_equal(before, after)


@dataclass(frozen=True)
class Verdict:
    """What `check_pass` found."""

    equivalent: bool
    dem_equal: bool
    witness: Optional[stim.Flow]

    @property
    def dem_blind(self) -> bool:
        """The DEM did not change, but the circuit did."""
        return self.dem_equal and not self.equivalent

    def __str__(self) -> str:
        if self.equivalent:
            return "equivalent"
        note = " (the DEM comparison misses this)" if self.dem_blind else ""
        found = f": {self.witness} holds in only one of them" if self.witness else ""
        return f"changed{note}{found}"


def check_pass(before, after, *, discard_final_state: bool = True) -> Verdict:
    """Decide whether a pass preserved its circuit, and whether a DEM check
    would have noticed if it did not. See `flows_equivalent` for the options."""
    equivalent, flow = flows_equivalent(
        before, after, discard_final_state=discard_final_state
    )
    return Verdict(equivalent, dem_equal(before, after), flow)


def witness() -> Tuple[stim.Circuit, stim.Circuit]:
    """The smallest pair this module was written for.

    One qubit, one detector, and a noise channel. The second circuit has an
    extra X between the two measurements, so its detector always fires. Both
    have the same DEM.
    """
    good = stim.Circuit("R 0\nX_ERROR(0.01) 0\nM 0\nM 0\nDETECTOR rec[-1] rec[-2]")
    bad = stim.Circuit("R 0\nX_ERROR(0.01) 0\nM 0\nX 0\nM 0\nDETECTOR rec[-1] rec[-2]")
    return good, bad


# --------------------------------------------------------------------------
# mutations: each one changes what a circuit computes, unless the spot it
# lands on happens to be dead. `oracle_power` uses the flows to tell which.
# --------------------------------------------------------------------------

_TWO_QUBIT = {"CX", "CY", "CZ", "XCZ", "YCZ", "SWAP", "ISWAP", "CXSWAP", "CZSWAP"}
_MEASURE = {"M", "MX", "MY", "MR", "MRX", "MRY"}


def _spots(c: stim.Circuit, names) -> list:
    return [i for i in range(len(c)) if isinstance(c[i], stim.CircuitInstruction) and c[i].name in names]


def _replace(c: stim.Circuit, i: int, new) -> stim.Circuit:
    out = stim.Circuit()
    for j in range(len(c)):
        if j == i:
            for ins in new:
                out.append(ins)
        else:
            out.append(c[j])
    return out


def _qubits(ins: stim.CircuitInstruction) -> list:
    return [t.value for t in ins.targets_copy() if t.is_qubit_target]


def mut_drop_two_qubit_gate(c: stim.Circuit, rng: random.Random) -> Optional[stim.Circuit]:
    """Remove one two-qubit gate instruction."""
    spots = _spots(c, _TWO_QUBIT)
    return _replace(c, rng.choice(spots), []) if spots else None


def mut_drop_control(c: stim.Circuit, rng: random.Random) -> Optional[stim.Circuit]:
    """Turn a CX into an X on its targets: the controls disappear."""
    spots = _spots(c, {"CX"})
    if not spots:
        return None
    i = rng.choice(spots)
    return _replace(c, i, [stim.CircuitInstruction("X", _qubits(c[i])[1::2])])


def mut_pauli_before_measurement(c: stim.Circuit, rng: random.Random) -> Optional[stim.Circuit]:
    """Insert a deterministic X right before a Z-basis measurement, flipping it.

    This is the mutation a DEM comparison cannot see.
    """
    spots = _spots(c, {"M", "MR"})
    if not spots:
        return None
    i = rng.choice(spots)
    q = rng.choice(_qubits(c[i]))
    return _replace(c, i, [stim.CircuitInstruction("X", [q]), c[i]])


def mut_swap_measurement_basis(c: stim.Circuit, rng: random.Random) -> Optional[stim.Circuit]:
    """Measure in X where the circuit measured in Z, or the other way round."""
    swap = {"M": "MX", "MX": "M", "MR": "MRX", "MRX": "MR"}
    spots = _spots(c, swap)
    if not spots:
        return None
    i = rng.choice(spots)
    ins = c[i]
    new = stim.CircuitInstruction(swap[ins.name], ins.targets_copy(), ins.gate_args_copy())
    return _replace(c, i, [new])


MUTATIONS = {
    "drop a two-qubit gate": mut_drop_two_qubit_gate,
    "drop the control of a CX": mut_drop_control,
    "X before a measurement": mut_pauli_before_measurement,
    "wrong measurement basis": mut_swap_measurement_basis,
}


@dataclass(frozen=True)
class Power:
    """How many changed circuits an oracle reported, for one mutation."""

    caught: int
    changed: int
    equivalent: int

    @property
    def rate(self) -> float:
        return self.caught / self.changed if self.changed else float("nan")


def oracle_power(
    oracle: Callable[[stim.Circuit, stim.Circuit], bool],
    circuits: Iterable,
    *,
    trials: int = 50,
    seed: int = 0,
    mutations=None,
    reference: Callable[[stim.Circuit, stim.Circuit], bool] = flows_differ,
) -> dict:
    """Measure how much of a known change an oracle can see.

    Each mutation is applied `trials` times to each circuit. Mutants that the
    reference says are unchanged (the mutation landed somewhere dead) are set
    aside, so the rate is over circuits that really changed.

    The reference is the flows check by default, so measuring `flows_differ`
    itself against it gives 100% by construction. Use this to measure other
    oracles: `dem_differs`, a detector sampler, or a project's own tests.

    Args:
        oracle: `oracle(before, after)` returns True when it reports a change.
            `dem_differs` and `flows_differ` are two such oracles.
        circuits: The circuits to mutate.
        trials: Mutants per mutation and circuit.
        seed: Seed for choosing where each mutation lands.
        mutations: Mapping of name to mutation, `MUTATIONS` by default.
        reference: Decides which mutants really changed the circuit.

    Returns:
        Mapping of mutation name to `Power`.
    """
    rng = random.Random(seed)
    mutations = MUTATIONS if mutations is None else mutations
    tally = {name: [0, 0, 0] for name in mutations}
    for circuit in circuits:
        base = to_stim(circuit).flattened()
        for name, mutate in mutations.items():
            for _ in range(trials):
                mutant = mutate(base, rng)
                if mutant is None:
                    break
                if not reference(base, mutant):
                    tally[name][2] += 1
                    continue
                tally[name][1] += 1
                if oracle(base, mutant):
                    tally[name][0] += 1
    return {name: Power(*counts) for name, counts in tally.items()}
