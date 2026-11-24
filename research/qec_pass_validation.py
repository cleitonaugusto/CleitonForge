"""Numbers for the note on validating QEC compiler passes.

For each circuit and mutation: how many mutants changed the circuit (by the
flows, confirmed independently by decoding the mutant's noiseless measurements
with the original's detectors and observables), and how many of those each
oracle reports.

    cforge-py/.venv/bin/python research/qec_pass_validation.py > research/qec_pass_validation.json
"""
import json
import random

import stim
from cleitonforge import qec

TRIALS = 50
NOISE = dict(after_clifford_depolarization=0.001, before_measure_flip_probability=0.001,
             after_reset_flip_probability=0.001)


def detectors_differ(a: stim.Circuit, b: stim.Circuit) -> bool:
    """Self-referential detector sampling: each circuit against its own reference."""
    if a.num_detectors != b.num_detectors:
        return True
    da = a.without_noise().compile_detector_sampler().sample(1)
    db = b.without_noise().compile_detector_sampler().sample(1)
    return bool((da != db).any())


def independently_changed(base: stim.Circuit, mutant: stim.Circuit) -> bool:
    """Decode the mutant's noiseless shots with the original's detectors and observables."""
    if mutant.num_measurements != base.num_measurements:
        return True
    conv = base.without_noise().compile_m2d_converter()
    shots = mutant.without_noise().compile_sampler().sample(64)
    det, obs = conv.convert(measurements=shots, separate_observables=True)
    return bool(det.any() or obs.any())


ORACLES = {"dem": qec.dem_differs, "detectors": detectors_differ, "flows": qec.flows_differ}


def circuits():
    for d in (3, 5):
        for task in ("surface_code:rotated_memory_z", "surface_code:rotated_memory_x",
                     "repetition_code:memory"):
            yield f"{task} d{d}", stim.Circuit.generated(task, distance=d, rounds=d, **NOISE)


def main():
    rows = []
    for name, c in circuits():
        base = c.flattened()
        for mname, mutate in qec.MUTATIONS.items():
            rng = random.Random(f"{name}/{mname}")
            row = {"circuit": name, "mutation": mname, "applied": 0, "changed": 0,
                   "independent_agrees": 0, **{o: 0 for o in ORACLES}}
            for _ in range(TRIALS):
                m = mutate(base, rng)
                if m is None:
                    break
                row["applied"] += 1
                changed = qec.flows_differ(base, m)
                row["independent_agrees"] += changed == independently_changed(base, m)
                if not changed:
                    continue
                row["changed"] += 1
                for oname, oracle in ORACLES.items():
                    row[oname] += oracle(base, m)
            rows.append(row)
    print(json.dumps({"trials": TRIALS, "noise": NOISE, "stim": stim.__version__, "rows": rows}, indent=1))


if __name__ == "__main__":
    main()
