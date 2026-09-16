"""Oracle power over QEC circuits, where there is no operator oracle.

`oracle_power_study.py` measures power on unitary circuits and decides
equivalence with the operator oracle, which is exact and therefore has power 1
by construction. A QEC circuit has resets and measurements, so it has no single
unitary and that construction does not carry over.

What replaces it here is the joint distribution of measurement outcomes. Two
circuits that agree on it are indistinguishable to anyone running them, which is
the honest definition of equivalence for this domain. It is estimated by
sampling rather than computed exactly, so it is declared as the reference oracle
rather than as ground truth, and every power below is relative to it. The
resolution of that reference is REFERENCE_RESOLUTION and is reported with the
table, because a power measured against a blunt reference is itself blunt.

Three oracles, in the order a QEC person would reach for them:

    measurements   the reference above
    detectors      distribution of detector outcomes, which are parities of
                   measurements, so strictly coarser
    dem            equality of the detector error model, which is what circuit
                   transformations in this field are usually validated with

Why the detector oracles come out at exactly zero, rather than merely low: a
detector is reported as a deviation from the value it takes when the circuit
runs without noise, and that reference is computed from the circuit itself. With
no noise, every detector is zero by construction, in both circuits, each against
its own reference. So a detector-based oracle is not weak here. It is empty, and
the emptiness is structural rather than a matter of sample size. That is the
regime in which a compiler pass has to be checked for semantic equivalence,
because the simplification has to happen before noise is added.

Two numbers per cell, and the distinction is the whole point:

    raw     detected / mutants applied
    power   detected / mutants that actually change the reference

The gap is the equivalent mutants. A pass that removes a genuinely dead gate
produces one on purpose, and it is kept in the mutation set as a negative
control: any oracle that flags it is reporting a fault that is not there.

    python3 qec_oracle_power.py [trials]
"""

from __future__ import annotations

import json
import os
import sys
from collections import defaultdict

import numpy as np
import stim

SHOTS = 40_000
REFERENCE_RESOLUTION = 8e-3   # two distributions closer than this read as equal
ORACLES = ["measurements", "detectors", "dem"]


# --------------------------------------------------------------------------
# circuits
# --------------------------------------------------------------------------
def repetition_round(n_data: int, rounds: int, rng=None) -> stim.Circuit:
    """A repetition-code memory experiment, written the long way.

    Deliberately not stim.Circuit.generated: the mutations below edit the gate
    list, and a hand-built circuit keeps that edit honest and readable.
    """
    data = list(range(0, 2 * n_data, 2))
    anc = list(range(1, 2 * n_data - 1, 2))
    c = stim.Circuit()
    c.append("R", data + anc)
    # Without this the circuit carries no information: every qubit starts in
    # |0>, nothing is noisy, and every measurement reads 0 whatever the mutation
    # did. The first version of this study reported every mutant as equivalent
    # for exactly that reason. A fixed random Pauli pattern gives the syndrome
    # something to be about.
    if rng is not None:
        flipped = [q for q in data if rng.random() < 0.4]
        if flipped:
            c.append("X", flipped)
    for r in range(rounds):
        for i, a in enumerate(anc):
            c.append("CX", [data[i], a])
            c.append("CX", [data[i + 1], a])
        c.append("M", anc)
        c.append("R", anc)
        for i in range(len(anc)):
            if r == 0:
                c.append("DETECTOR", [stim.target_rec(-len(anc) + i)])
            else:
                c.append("DETECTOR", [stim.target_rec(-len(anc) + i),
                                      stim.target_rec(-2 * len(anc) + i)])
    c.append("M", data)
    return c


# --------------------------------------------------------------------------
# mutations, chosen to look like things a compiler pass does
# --------------------------------------------------------------------------
def _instructions(c: stim.Circuit) -> list:
    return [c[i] for i in range(len(c))]


def _rebuild(instrs) -> stim.Circuit:
    out = stim.Circuit()
    for ins in instrs:
        out.append(ins)
    return out


def mut_drop_dead_gate(c, rng):
    """Insert then drop an H immediately before a reset: an equivalent mutant.

    This is the negative control. The gate it removes cannot be observed, so an
    oracle that reports a difference here is wrong.
    """
    instrs = _instructions(c)
    spots = [i for i, ins in enumerate(instrs) if ins.name == "R"]
    if not spots:
        return None
    i = int(rng.choice(spots))
    targets = [t.value for t in instrs[i].targets_copy()]
    with_dead = instrs[:i] + [stim.CircuitInstruction("H", targets)] + instrs[i:]
    return _rebuild(with_dead), _rebuild(instrs)


def mut_drop_live_gate(c, rng):
    """Remove a two-qubit gate. Almost always changes the circuit."""
    instrs = _instructions(c)
    spots = [i for i, ins in enumerate(instrs) if ins.name == "CX"]
    if not spots:
        return None
    i = int(rng.choice(spots))
    return _rebuild(instrs), _rebuild(instrs[:i] + instrs[i + 1:])


def mut_drop_control(c, rng):
    """Turn a CX into an X on its target: the control silently disappears."""
    instrs = _instructions(c)
    spots = [i for i, ins in enumerate(instrs) if ins.name == "CX"]
    if not spots:
        return None
    i = int(rng.choice(spots))
    tg = [t.value for t in instrs[i].targets_copy()]
    replaced = stim.CircuitInstruction("X", tg[1::2])
    return _rebuild(instrs), _rebuild(instrs[:i] + [replaced] + instrs[i + 1:])


def mut_insert_deterministic_pauli(c, rng):
    """Insert an X between a measurement and the next reset on that qubit."""
    instrs = _instructions(c)
    spots = [i for i, ins in enumerate(instrs) if ins.name == "M"]
    if not spots:
        return None
    i = int(rng.choice(spots))
    tg = [t.value for t in instrs[i].targets_copy()]
    inserted = stim.CircuitInstruction("X", [tg[0]])
    return _rebuild(instrs), _rebuild(instrs[:i + 1] + [inserted] + instrs[i + 1:])


MUTATIONS = {
    "drop dead gate (equivalent, control)": mut_drop_dead_gate,
    "drop a CX": mut_drop_live_gate,
    "drop the control of a CX": mut_drop_control,
    # Lands between a measurement and the reset of that same ancilla, so the
    # Pauli it inserts is genuinely dead. It came out equivalent on every trial
    # and is kept as a second negative control rather than relabelled.
    "insert a deterministic Pauli (equivalent, control)": mut_insert_deterministic_pauli,
}


# --------------------------------------------------------------------------
# oracles
# --------------------------------------------------------------------------
def _dist(sample: np.ndarray) -> np.ndarray:
    return sample.mean(axis=0)


def oracle_measurements(a: stim.Circuit, b: stim.Circuit, rng) -> bool:
    """The reference. True means the two circuits look different.

    These circuits carry no noise, so every measurement is deterministic and one
    shot settles it exactly. That is worth more than a large sample: the power
    figures below are relative to an exact reference, not an estimated one.
    """
    sa = a.compile_sampler().sample(1)
    sb = b.compile_sampler().sample(1)
    if sa.shape[1] != sb.shape[1]:
        return True
    return bool((sa != sb).any())


def oracle_detectors(a: stim.Circuit, b: stim.Circuit, rng) -> bool:
    seed = int(rng.integers(1 << 30))
    da = a.compile_detector_sampler().sample(1)
    db = b.compile_detector_sampler().sample(1)
    if da.shape[1] != db.shape[1]:
        return True
    return bool((da != db).any())


def oracle_dem(a: stim.Circuit, b: stim.Circuit, rng) -> bool:
    try:
        return str(a.detector_error_model(decompose_errors=False)) != str(
            b.detector_error_model(decompose_errors=False))
    except Exception:
        # A mutant whose DEM will not build is a difference the oracle can see.
        return True


ORACLE_FNS = {
    "measurements": oracle_measurements,
    "detectors": oracle_detectors,
    "dem": oracle_dem,
}


# --------------------------------------------------------------------------
# study
# --------------------------------------------------------------------------
def run_cell(rng, n_data, rounds, mut_fn, trials):
    applied = equivalent = false_pos = 0
    detected = defaultdict(int)

    for _ in range(trials):
        base = repetition_round(n_data, rounds, rng)
        made = mut_fn(base, rng)
        if made is None:
            continue
        good, bad = made
        applied += 1

        if not oracle_measurements(good, bad, rng):
            equivalent += 1
            for name in ("detectors", "dem"):
                if ORACLE_FNS[name](good, bad, rng):
                    false_pos += 1
                    break
            continue

        detected["measurements"] += 1        # true by construction, it decided
        for name in ("detectors", "dem"):
            if ORACLE_FNS[name](good, bad, rng):
                detected[name] += 1

    return applied, equivalent, false_pos, dict(detected)


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rng = np.random.default_rng(20260916)
    rows = []

    for n_data in (3, 5):
        for rounds in (2, 4):
            for mut_name, mut_fn in MUTATIONS.items():
                applied, equiv, fp, det = run_cell(
                    rng, n_data, rounds, mut_fn, trials)
                changed = applied - equiv
                rows.append({
                    "n_data": n_data, "rounds": rounds, "mutation": mut_name,
                    "applied": applied, "equivalent": equiv, "false_positives": fp,
                    "oracles": {
                        o: {
                            "detected": det.get(o, 0),
                            "raw": det.get(o, 0) / applied if applied else None,
                            "power": det.get(o, 0) / changed if changed else None,
                        } for o in ORACLES
                    },
                })

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "qec_oracle_power.json")
    with open(out, "w") as fh:
        json.dump({"shots": SHOTS, "resolution": REFERENCE_RESOLUTION,
                   "trials_per_cell": trials, "rows": rows}, fh, indent=2)

    print(f"{'mutation':<38} {'applied':>7} {'equiv':>6} "
          + " ".join(f"{o:>13}" for o in ORACLES))
    print("-" * 96)
    for mut in MUTATIONS:
        sel = [r for r in rows if r["mutation"] == mut]
        ap = sum(r["applied"] for r in sel)
        eq = sum(r["equivalent"] for r in sel)
        ch = ap - eq
        cells = []
        for o in ORACLES:
            d = sum(r["oracles"][o]["detected"] for r in sel)
            cells.append(f"{d}/{ch}" if ch else "n/a")
        print(f"{mut:<38} {ap:>7} {eq:>6} " + " ".join(f"{c:>13}" for c in cells))
    print("-" * 96)
    print(f"  power = detected / mutants that changed the reference")
    print(f"  reference = exact measurement record (noiseless, deterministic)")
    print(f"  wrote {os.path.basename(out)}")


if __name__ == "__main__":
    main()
