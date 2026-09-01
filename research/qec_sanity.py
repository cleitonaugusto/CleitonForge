"""Known-answer battery for a QEC reference oracle. Run this before measuring.

Four reference-oracle designs were tried on 2026-09-01 and all four were wrong,
each printing a plausible number. The common cause was that every one of them
was *self-referential*: it let each circuit define its own ground truth. stim
derives the detector reference sample from whichever circuit you sample, so a
deterministic fault shifts the reference along with the outcome and vanishes.
A differential test has to judge the mutant against the original's reference,
which is what compile_m2d_converter does.

This file asserts what the answers must be, on cases where the answer is known
independently of any measurement. An oracle that fails any of these measures
nothing, whatever it prints.

    .venv-qec/bin/python qec_sanity.py
"""
import sys

import numpy as np
import stim

SHOTS = 20000
SEED = 20260901

DISTANCE, ROUNDS, NOISE = 3, 3, 0.002


def base_circuit(noise=NOISE):
    return stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=DISTANCE, rounds=ROUNDS,
        after_clifford_depolarization=noise,
        after_reset_flip_probability=noise,
        before_measure_flip_probability=noise,
        before_round_data_depolarization=noise,
    )


def insert(circuit, gate, qubit, at):
    """Insert a single-qubit gate before flattened instruction index `at`."""
    out = stim.Circuit()
    for j, inst in enumerate(circuit.flattened()):
        if j == at:
            out.append(gate, [qubit])
        out.append(inst)
    return out


# --------------------------------------------------------------- the oracle
def cross_referential(original, mutant, shots=SHOTS, seed=SEED):
    """Detection-event rates of the mutant, judged by the original's definitions.

    This is the whole fix. `compile_m2d_converter` takes the original circuit's
    detector definitions and reference sample and applies them to somebody
    else's measurement records, so a deterministic fault can no longer hide by
    dragging the reference along with it.

    Returns (detector_rates, observable_rate) or None if the shapes do not line
    up, which is itself a difference and a loud one.
    """
    conv = original.compile_m2d_converter()
    records = mutant.compile_sampler(seed=seed).sample(shots)
    if records.shape[1] != original.num_measurements:
        return None
    det, obs = conv.convert(measurements=records, separate_observables=True)
    return det.mean(axis=0), obs.mean(axis=0)


def differs(a, b, shots=SHOTS, z=5.0):
    """Two-proportion test, per detector and on the observable."""
    if a is None or b is None:
        return True
    for pa, pb in ((a[0], b[0]), (a[1], b[1])):
        pool = (pa + pb) / 2
        se = np.sqrt(np.maximum(pool * (1 - pool), 1e-12) * 2 / shots)
        if np.any(np.abs(pa - pb) / se > z):
            return True
    return False


# --------------------------------------------------------------- known answers
def find_sites(circuit):
    """A reset site and a data qubit, both needed for the known-answer cases."""
    flat = list(circuit.flattened())
    reset_site = None
    for i, inst in enumerate(flat):
        if inst.name in ("R", "RX", "MR"):
            qs = [t.value for t in inst.targets_copy() if t.is_qubit_target]
            if qs:
                reset_site = (i + 1, qs[0])   # immediately after the reset
                break

    # data qubits are the ones measured only at the very end
    last_m = max(i for i, inst in enumerate(flat) if inst.name == "M")
    data = sorted({t.value for t in flat[last_m].targets_copy()
                   if t.is_qubit_target})
    return reset_site, data, last_m


def main():
    base = base_circuit()
    noiseless = base_circuit(noise=0.0)
    (reset_at, reset_q), data, last_m = find_sites(base)
    mid = len(list(base.flattened())) // 2

    print(f"rotated surface code d={DISTANCE}, {ROUNDS} rounds, {SHOTS} shots")
    print(f"reset site: instruction {reset_at}, qubit {reset_q}")
    print(f"data qubits: {data}\n")

    ref = cross_referential(base, base)

    cases = []

    # --- must be EQUIVALENT ------------------------------------------------
    cases.append((
        "Z right after a reset",
        "Z|0> = |0>, so this gate is the identity",
        insert(base, "Z", reset_q, reset_at), False))

    cases.append((
        "Z twice on the same qubit",
        "ZZ = I",
        insert(insert(base, "Z", data[0], mid), "Z", data[0], mid), False))

    cases.append((
        "S then S_DAG",
        "S S+ = I",
        insert(insert(base, "S_DAG", data[0], mid), "S", data[0], mid), False))

    cases.append((
        "the circuit against itself",
        "identical input must not raise an alarm",
        base.copy(), False))

    # --- must DIFFER -------------------------------------------------------
    cases.append((
        "X on a data qubit before the final measurement",
        "flips that measurement, so the Z detectors touching it must fire",
        insert(base, "X", data[0], last_m), True))

    cases.append((
        "X on a data qubit mid-circuit",
        "a bit flip is exactly what this code is built to detect",
        insert(base, "X", data[0], mid), True))

    cases.append((
        "Z on a data qubit mid-circuit",
        "anticommutes with the X stabilisers touching that qubit",
        insert(base, "Z", data[0], mid), True))

    # --- run ---------------------------------------------------------------
    width = max(len(name) for name, _, _, _ in cases) + 2
    passed = failed = 0
    print(f"{'case':<{width}}{'expected':>12}{'got':>10}   verdict")
    print("-" * (width + 40))
    for name, why, mutant, should_differ in cases:
        got = differs(ref, cross_referential(base, mutant))
        ok = (got == should_differ)
        passed, failed = passed + ok, failed + (not ok)
        print(f"{name:<{width}}{'differs' if should_differ else 'equivalent':>12}"
              f"{'differs' if got else 'equivalent':>10}   "
              f"{'ok' if ok else 'FAIL'}")
        if not ok:
            print(f"{'':<{width}}  why: {why}")

    print(f"\n{passed} passed, {failed} failed")
    if failed:
        print("\nThis oracle is not usable. Fix it before measuring anything.")
        return 1
    print("\nThe oracle answers every known case correctly. It may now be used")
    print("to measure power on cases whose answer is not known.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
