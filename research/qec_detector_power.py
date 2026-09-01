"""WORK IN PROGRESS -- the numbers this prints are NOT trustworthy yet.

Four reference-oracle designs were tried on 2026-09-01 and all four were wrong.
See CForge-Vault/Pesquisa/QEC-Porte-Tentativa-01-09.md. In short: comparing raw
DEM text reports bookkeeping as physics; comparing error mechanisms alone calls
every deterministic fault equivalent; and stim's detector sampler reports events
relative to each circuit's own reference sample, so it absorbs deterministic
changes entirely. Do not cite this file until a known-answer sanity battery
passes -- Z right after a reset must come out equivalent, an X on the
observable's support must flip it.

Power of the detector-sample oracle on a real error-correction circuit.

The Clifford blindness measured in oracle_power_study.py was about outcome
counts from |0...0> in random circuits. Error-correction tooling is a different
setting -- measurements, resets, detectors, noise -- so whether the blindness
carries over is an empirical question, not a corollary. This measures it.

Target: stim's rotated surface code memory experiment, which is the workload a
QEC compiler test suite actually runs.

Faults injected, one gate at a time, at every legal site:

  insert S / insert Z / insert H   a spurious gate a pass might leave behind
  swap CX control and target       the classic qubit-ordering bug
  drop CX / drop H                 an over-eager cancellation pass

Two oracles, the same design as the main study:

  strong   the detector error model, compared after canonicalisation. This is
           the reference: mutants it cannot separate are equivalent, and no
           correct oracle can flag them.
  weak     detector-sample frequencies over N shots, compared per detector with
           a two-proportion test. This is what a test suite does.

Power is detected / non-equivalent, and the weak oracle's false-positive rate
is measured on the equivalent mutants, which are a free negative control.

    .venv-qec/bin/python qec_detector_power.py [shots]
"""
import sys
from collections import defaultdict

import numpy as np
import stim

SHOTS = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
SEED = 20260901

DISTANCE, ROUNDS = 3, 3
NOISE = 0.002


def base_circuit():
    return stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=DISTANCE, rounds=ROUNDS,
        after_clifford_depolarization=NOISE,
        after_reset_flip_probability=NOISE,
        before_measure_flip_probability=NOISE,
        before_round_data_depolarization=NOISE,
    )


# ------------------------------------------------------------------ mutations
def mutants(circuit):
    """Yield (fault_class, site, mutated_circuit) for every legal injection."""
    flat = list(circuit.flattened())

    for i, inst in enumerate(flat):
        name = inst.name

        # drop a gate
        if name in ("CX", "H"):
            out = stim.Circuit()
            for j, other in enumerate(flat):
                if j != i:
                    out.append(other)
            yield (f"drop {name}", i, out)

        # swap the two qubits of a CX
        if name == "CX":
            targets = inst.targets_copy()
            if len(targets) >= 2:
                swapped = list(targets)
                swapped[0], swapped[1] = swapped[1], swapped[0]
                out = stim.Circuit()
                for j, other in enumerate(flat):
                    if j == i:
                        out.append("CX", swapped)
                    else:
                        out.append(other)
                yield ("swap CX qubits", i, out)

    # Insert a spurious single-qubit gate. Enumerate over every qubit at every
    # TICK boundary, not over instruction targets: taking the first target of
    # each instruction looks like many sites but lands on a handful of qubits,
    # which silently turns a survey into a spot check.
    ticks = [i for i, inst in enumerate(flat) if inst.name == "TICK"]
    qubits = sorted(circuit.get_final_qubit_coordinates())
    for i in ticks:
        for q in qubits:
            for g in ("S", "Z", "H"):
                out = stim.Circuit()
                for j, other in enumerate(flat):
                    if j == i:
                        out.append(g, [q])
                    out.append(other)
                yield (f"insert {g}", (i, q), out)


# ------------------------------------------------------------------ oracles
def dem_key(circuit):
    """Semantic form of the detector error model, or a marker if it will not build.

    Only the error mechanisms are kept: each `error(p) D.. L..` line, with its
    targets sorted, and the whole multiset sorted. Everything else in a DEM is
    bookkeeping -- `detector(x,y,t)` coordinate declarations and
    `shift_detectors` folding -- and inserting an instruction perturbs that
    bookkeeping without changing any physics. Comparing the raw text therefore
    reports a difference for a gate that is provably the identity, which
    inflates the non-equivalent denominator and makes every weaker oracle look
    blind. Detector indices are comparable because the mutations here do not
    add or remove detectors; that is asserted by the caller.
    """
    try:
        dem = circuit.detector_error_model(
            decompose_errors=False, allow_gauge_detectors=True).flattened()
    except Exception:
        return "UNBUILDABLE"
    mechanisms = []
    for line in str(dem).splitlines():
        line = line.strip()
        if not line.startswith("error("):
            continue
        head, _, targets = line.partition(") ")
        mechanisms.append(head + ") " + " ".join(sorted(targets.split())))

    # The error mechanisms are only half of it. A deterministic fault -- a
    # stray Z, say -- adds no noise term at all; it moves the noiseless
    # reference the detectors are defined against. Comparing mechanisms alone
    # calls every such fault equivalent, which is as wrong in the other
    # direction as comparing raw text. So the noiseless signature goes in too.
    try:
        clean = circuit.without_noise()
        sample = clean.compile_detector_sampler(seed=1).sample(
            1, separate_observables=True)
        det, obs = sample
        signature = ("".join("1" if b else "0" for b in det[0]) + "|" +
                     "".join("1" if b else "0" for b in obs[0]))
    except Exception:
        signature = "NOISELESS-UNBUILDABLE"

    return signature + "\n" + "\n".join(sorted(mechanisms))


def sample_detectors(circuit, shots, seed):
    sampler = circuit.compile_detector_sampler(seed=seed)
    return sampler.sample(shots, bit_packed=False)


def weak_oracle(freq_a, freq_b, shots, z=4.0):
    """Two-proportion z test per detector, Bonferroni-style cut at z=4.

    Any detector whose firing rate differs by more than the cut counts as a
    detection. z=4 keeps the per-circuit false-positive rate low across the
    ~24 detectors; the actual rate is measured, not assumed.
    """
    n = shots
    pa, pb = freq_a, freq_b
    pool = (pa + pb) / 2
    se = np.sqrt(np.maximum(pool * (1 - pool), 1e-12) * 2 / n)
    return bool(np.any(np.abs(pa - pb) / se > z))


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (p, max(0.0, c - h), min(1.0, c + h))


# ------------------------------------------------------------------ main
def main():
    base = base_circuit()
    print(f"rotated surface code d={DISTANCE}, {ROUNDS} rounds, noise {NOISE}")
    print(f"{base.num_detectors} detectors, {base.num_qubits} qubits, "
          f"{SHOTS} shots per circuit\n")

    base_dem = dem_key(base)
    base_freq = sample_detectors(base, SHOTS, SEED).mean(axis=0)

    stats = defaultdict(lambda: {"gen": 0, "equiv": 0, "det": 0, "fp": 0,
                                 "unbuildable": 0})
    seed = SEED

    for fault, site, mut in mutants(base):
        seed += 1
        s = stats[fault]
        s["gen"] += 1

        mdem = dem_key(mut)
        if mdem == "UNBUILDABLE":
            # the fault makes the circuit ill-formed; a compiler bug that does
            # this is caught by construction, so it is not an oracle question
            s["unbuildable"] += 1
            continue

        equivalent = (mdem == base_dem)
        try:
            freq = sample_detectors(mut, SHOTS, seed).mean(axis=0)
        except Exception:
            s["unbuildable"] += 1
            continue

        if freq.shape != base_freq.shape:
            # different detector count is itself a detection, and a loud one
            if equivalent:
                s["equiv"] += 1
            else:
                s["det"] += 1
            continue

        flagged = weak_oracle(base_freq, freq, SHOTS)
        if equivalent:
            s["equiv"] += 1
            if flagged:
                s["fp"] += 1
        elif flagged:
            s["det"] += 1

    # -------------------------------------------------------------- report
    print(f"{'fault class':<20}{'sites':>7}{'ill-formed':>12}"
          f"{'equivalent':>12}{'power':>18}")
    print("-" * 69)
    tot = defaultdict(int)
    for fault in sorted(stats):
        s = stats[fault]
        live = s["gen"] - s["equiv"] - s["unbuildable"]
        p, lo, hi = wilson(s["det"], live)
        for k in s:
            tot[k] += s[k]
        tot["live"] += live
        tot["detsum"] += s["det"]
        cell = (f"{p*100:>8.1f}%  [{lo*100:.0f},{hi*100:.0f}]"
                if live else "        --        ")
        print(f"{fault:<20}{s['gen']:>7}{s['unbuildable']:>12}"
              f"{s['equiv']:>12}{cell:>18}")

    p, lo, hi = wilson(tot["detsum"], tot["live"])
    print("-" * 69)
    print(f"{'ALL':<20}{tot['gen']:>7}{tot['unbuildable']:>12}"
          f"{tot['equiv']:>12}{p*100:>8.1f}%  [{lo*100:.0f},{hi*100:.0f}]")

    if tot["equiv"]:
        fp, flo, fhi = wilson(tot["fp"], tot["equiv"])
        print(f"\nweak-oracle false positives on the {tot['equiv']} equivalent "
              f"mutants: {tot['fp']} = {fp*100:.1f}% [{flo*100:.1f}, {fhi*100:.1f}]")
    print("\npower = detected / (sites - equivalent - ill-formed)")
    print("equivalent = the detector error model cannot separate it from the "
          "original,\n             so no correct oracle can flag it")


if __name__ == "__main__":
    main()
