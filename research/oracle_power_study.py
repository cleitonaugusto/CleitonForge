"""Oracle power across circuit families, widths, depths and fault classes.

The question is not "how many bugs did we find". It is "what fraction of
injected faults could this oracle have detected at all", measured before any
claim about a target.

Two numbers per cell, and the distinction is the point:

  raw     detected / mutants generated
  power   detected / mutants that actually change the circuit

The gap between them is the equivalent mutants: faults that produce an
operator identical up to global phase. No correct oracle can flag those, and
counting them as misses understates every oracle. Equivalence here is decided
by the operator oracle, which therefore has power 1 by construction -- it is
the reference, not a competitor. What the study measures is how far the
cheaper oracles fall below it.

    python3 oracle_power_study.py [trials]

Writes oracle_power_study.json next to itself.
"""
import json
import sys
import time
from collections import defaultdict

import numpy as np

import qcore as q

TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 400
SEED = 20260901

FAMILIES = {
    "blind {H,X,Z,Ry,Rz,CX,CZ}": ["H", "X", "Z", "Ry", "Rz", "CX", "CZ"],
    "Clifford {H,S,CX,CZ}": ["H", "S", "CX", "CZ"],
    "Clifford+T": ["H", "S", "T", "CX", "CZ"],
    "variational {Rx,Ry,Rz,CX}": ["Rx", "Ry", "Rz", "CX"],
}
WIDTHS = [3, 4, 5]
DEPTHS = [10, 20, 40]
ORACLES = ["sampled counts (1024 shots)", "exact counts", "state", "operator"]


def run_cell(rng, gates, n, depth, mut_name, mut_fn, trials):
    applied = 0
    equivalent = 0
    false_pos = 0
    detected = defaultdict(int)

    for _ in range(trials):
        ops = q.random_circuit(rng, n, depth, gates)
        mutant = mut_fn(ops, rng)
        if mutant is None:
            continue
        applied += 1

        Ua = q.simulate(ops, n)
        Ub = q.simulate(mutant, n)

        if not q.oracle_operator(Ua, Ub):
            equivalent += 1
            # An equivalent mutant is a free negative control: the operators
            # agree, so any alarm here is a false positive. Only the sampled
            # oracle can produce one, and it does -- which is why its raw
            # detection rate is not a power.
            if q.oracle_sampled_counts(Ua, Ub, rng):
                false_pos += 1
            continue

        if q.oracle_sampled_counts(Ua, Ub, rng):
            detected["sampled counts (1024 shots)"] += 1
        if q.oracle_exact_counts(Ua, Ub):
            detected["exact counts"] += 1
        if q.oracle_state(Ua, Ub):
            detected["state"] += 1
        detected["operator"] += 1          # true by construction, see docstring

    return applied, equivalent, false_pos, dict(detected)


def main():
    rng = np.random.default_rng(SEED)
    t0 = time.time()
    rows = []

    total = len(FAMILIES) * len(WIDTHS) * len(DEPTHS) * len(q.MUTATIONS)
    done = 0

    for fam, gates in FAMILIES.items():
        for n in WIDTHS:
            for depth in DEPTHS:
                for mut_name, mut_fn in q.MUTATIONS.items():
                    applied, equiv, fp, det = run_cell(
                        rng, gates, n, depth, mut_name, mut_fn, TRIALS)
                    done += 1
                    if applied == 0:
                        continue
                    live = applied - equiv
                    row = {
                        "family": fam, "width": n, "depth": depth,
                        "mutation": mut_name,
                        "generated": applied, "equivalent": equiv,
                        "non_equivalent": live,
                        "sampled_false_positives": fp,
                        "oracles": {},
                    }
                    for o in ORACLES:
                        k = det.get(o, 0)
                        p, lo, hi = q.wilson(k, live)
                        row["oracles"][o] = {
                            "detected": k, "power": p, "ci95": [lo, hi]}
                    rows.append(row)
                    print(f"\r  {done}/{total} cells  "
                          f"({time.time()-t0:.0f}s)", end="", flush=True)

    print()
    out = {"trials_per_cell": TRIALS, "seed": SEED,
           "seconds": round(time.time() - t0, 1), "rows": rows}
    with open("oracle_power_study.json", "w") as f:
        json.dump(out, f, indent=1)

    report(rows)


def agg(rows, key):
    """Pool cells by one factor and recompute Wilson on the pooled counts."""
    buckets = defaultdict(lambda: defaultdict(int))
    for r in rows:
        b = buckets[r[key]]
        b["live"] += r["non_equivalent"]
        b["gen"] += r["generated"]
        b["equiv"] += r["equivalent"]
        b["fp"] += r["sampled_false_positives"]
        for o in ORACLES:
            b[o] += r["oracles"][o]["detected"]
    return buckets


def table(title, buckets, label_w=28):
    print(f"\n{title}")
    head = f"{'':<{label_w}}{'equiv':>7}" + "".join(
        f"{o.split(' (')[0]:>16}" for o in ORACLES[:3])
    print(head)
    print("-" * len(head))
    for name, b in buckets.items():
        live, gen = b["live"], b["gen"]
        eq = b["equiv"] / gen * 100 if gen else 0
        cells = ""
        for o in ORACLES[:3]:
            p, lo, hi = q.wilson(b[o], live)
            cells += f"{p*100:>9.1f}% ±{(hi-lo)/2*100:>4.1f}"
        print(f"{str(name):<{label_w}}{eq:>6.1f}%{cells}")


def report(rows):
    tot_gen = sum(r["generated"] for r in rows)
    tot_eq = sum(r["equivalent"] for r in rows)
    tot_live = tot_gen - tot_eq
    print(f"\n{tot_gen} mutants generated, {tot_eq} equivalent "
          f"({tot_eq/tot_gen*100:.1f}%), {tot_live} non-equivalent.")
    tot_fp = sum(r["sampled_false_positives"] for r in rows)
    p, lo, hi = q.wilson(tot_fp, tot_eq)
    print("Power below is against non-equivalent mutants only. "
          "The operator oracle is the reference (power 1 by construction) "
          "and is omitted from the tables.")
    print(f"\nSampled-counts false positives, measured on the {tot_eq} "
          f"equivalent mutants (operators identical, so every alarm is wrong):"
          f"\n  {tot_fp}/{tot_eq} = {p*100:.1f}%  "
          f"[95% CI {lo*100:.1f}%, {hi*100:.1f}%]"
          f"\n  Its detection rate is inflated by roughly this much and is "
          f"therefore not a power.")

    for key, title in [("mutation", "By fault class"),
                       ("family", "By circuit family"),
                       ("width", "By width (qubits)"),
                       ("depth", "By depth (gates)")]:
        table(title, agg(rows, key))


if __name__ == "__main__":
    main()
