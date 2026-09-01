"""Power of the detector-sample oracle on a real error-correction circuit.

Run qec_sanity.py first. It asserts the known answers this file's reference
oracle has to get right, and four earlier designs did not.

The Clifford blindness in oracle_power_study.py was about outcome counts from
|0...0> in random circuits. Error correction is a different setting, so whether
it carries over is an empirical question rather than a corollary.

Three oracles, and the contrast between the first two and the third is the
result:

  self-referential detectors    each circuit sampled with its own detector
                                sampler, detection-event rates compared. This
                                is what a QEC test suite does.
  self-referential + logical    the same, plus the observable's flip rate.
  cross-referential             the mutant's measurement records judged by the
                                ORIGINAL circuit's detector definitions and
                                reference sample, via compile_m2d_converter.

A self-referential oracle derives its ground truth from the circuit under test.
stim builds the detector reference sample from whichever circuit you hand it,
so a fault that shifts the outcome shifts the reference along with it and
leaves no trace. The cross-referential oracle is the differential one, and is
therefore the reference: equivalence is decided by it.

    .venv-qec/bin/python qec_detector_power.py [shots]
"""
import sys
from collections import defaultdict

import numpy as np
import stim

from qec_sanity import base_circuit, cross_referential, differs

SHOTS = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
SEED = 20260901


def self_referential(circuit, shots, seed):
    """Detection-event and observable rates from the circuit's own sampler."""
    det, obs = circuit.compile_detector_sampler(seed=seed).sample(
        shots, separate_observables=True)
    return det.mean(axis=0), obs.mean(axis=0)


def rates_differ(a, b, shots, z=5.0, with_observable=True):
    if a is None or b is None:
        return True
    pairs = [(a[0], b[0])] + ([(a[1], b[1])] if with_observable else [])
    for pa, pb in pairs:
        if pa.shape != pb.shape:
            return True
        pool = (pa + pb) / 2
        se = np.sqrt(np.maximum(pool * (1 - pool), 1e-12) * 2 / shots)
        if np.any(np.abs(pa - pb) / se > z):
            return True
    return False


def mutants(circuit):
    flat = list(circuit.flattened())
    ticks = [i for i, inst in enumerate(flat) if inst.name == "TICK"]
    qubits = sorted(circuit.get_final_qubit_coordinates())

    for i, inst in enumerate(flat):
        if inst.name in ("CX", "H"):
            out = stim.Circuit()
            for j, other in enumerate(flat):
                if j != i:
                    out.append(other)
            yield (f"drop {inst.name}", out)

        if inst.name == "CX":
            t = list(inst.targets_copy())
            if len(t) >= 2:
                t[0], t[1] = t[1], t[0]
                out = stim.Circuit()
                for j, other in enumerate(flat):
                    out.append(stim.CircuitInstruction("CX", t)
                               if j == i else other)
                yield ("swap CX qubits", out)

    # Every qubit at every TICK. Enumerating over instruction targets instead
    # looks like many sites but lands on a handful of qubits.
    for i in ticks:
        for q in qubits:
            for g in ("X", "Z", "S", "H"):
                out = stim.Circuit()
                for j, other in enumerate(flat):
                    if j == i:
                        out.append(g, [q])
                    out.append(other)
                yield (f"insert {g}", out)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (p, max(0.0, c - h), min(1.0, c + h))


def main():
    base = base_circuit()
    print(f"{base.num_detectors} detectors, {base.num_qubits} qubits, "
          f"{SHOTS} shots per circuit")
    print("power is against mutants the cross-referential oracle separates\n")

    ref_cross = cross_referential(base, base, SHOTS, SEED)
    ref_self = self_referential(base, SHOTS, SEED)

    stats = defaultdict(lambda: defaultdict(int))
    seed = SEED

    for fault, mut in mutants(base):
        seed += 1
        s = stats[fault]
        s["gen"] += 1

        try:
            cross = cross_referential(base, mut, SHOTS, seed)
        except Exception:
            s["ill"] += 1
            continue

        if not differs(ref_cross, cross, SHOTS):
            s["equiv"] += 1
            # equivalent mutants are a free negative control for the weak ones
            try:
                m_self = self_referential(mut, SHOTS, seed)
                if rates_differ(ref_self, m_self, SHOTS, with_observable=False):
                    s["fp_det"] += 1
                if rates_differ(ref_self, m_self, SHOTS):
                    s["fp_log"] += 1
            except Exception:
                pass
            continue

        s["live"] += 1
        try:
            m_self = self_referential(mut, SHOTS, seed)
        except Exception:
            continue
        if rates_differ(ref_self, m_self, SHOTS, with_observable=False):
            s["det"] += 1
        if rates_differ(ref_self, m_self, SHOTS):
            s["log"] += 1

    hdr = (f"{'fault class':<20}{'sites':>7}{'ill':>5}{'equiv':>7}{'live':>6}"
           f"{'self-ref detectors':>22}{'+ observable':>18}")
    print(hdr)
    print("-" * len(hdr))
    tot = defaultdict(int)
    for fault in sorted(stats):
        s = stats[fault]
        for k, v in s.items():
            tot[k] += v
        live = s["live"]
        if not live:
            print(f"{fault:<20}{s['gen']:>7}{s['ill']:>5}{s['equiv']:>7}"
                  f"{live:>6}{'--':>22}{'--':>18}")
            continue
        pd, dlo, dhi = wilson(s["det"], live)
        pl, llo, lhi = wilson(s["log"], live)
        print(f"{fault:<20}{s['gen']:>7}{s['ill']:>5}{s['equiv']:>7}{live:>6}"
              f"{pd*100:>13.1f}% [{dlo*100:>3.0f},{dhi*100:>3.0f}]"
              f"{pl*100:>9.1f}% [{llo*100:>3.0f},{lhi*100:>3.0f}]")

    print("-" * len(hdr))
    pd, dlo, dhi = wilson(tot["det"], tot["live"])
    pl, llo, lhi = wilson(tot["log"], tot["live"])
    print(f"{'ALL':<20}{tot['gen']:>7}{tot['ill']:>5}{tot['equiv']:>7}"
          f"{tot['live']:>6}{pd*100:>13.1f}% [{dlo*100:>3.0f},{dhi*100:>3.0f}]"
          f"{pl*100:>9.1f}% [{llo*100:>3.0f},{lhi*100:>3.0f}]")

    if tot["equiv"]:
        f1, a, b = wilson(tot["fp_det"], tot["equiv"])
        f2, c, d = wilson(tot["fp_log"], tot["equiv"])
        print(f"\nfalse positives on the {tot['equiv']} equivalent mutants: "
              f"detectors {f1*100:.1f}% [{a*100:.1f},{b*100:.1f}], "
              f"with observable {f2*100:.1f}% [{c*100:.1f},{d*100:.1f}]")


if __name__ == "__main__":
    main()
