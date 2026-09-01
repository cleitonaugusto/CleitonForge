"""Power hierarchy: counts < state < operator, for the Rz convention fault.

This reproduces the table published in README.md. It is deliberately narrow --
one fault class, one width, one depth -- because it exists to make the blindness
theorem visible in the smallest setting where it is true. The factorial version
lives in oracle_power_study.py.

The fault is a *convention* bug: every Rz in the circuit is built with the
opposite sign. That is not the same as flipping one Rz, and the difference is
the whole result. Inverting every Rz conjugates a circuit made of otherwise
real gates, and |z|^2 = |conj(z)|^2, so no measurement probability moves.
Inverting a single Rz conjugates nothing and is ordinarily visible.

Power is measured against non-equivalent mutants: those the operator oracle
can tell apart. Counting equivalent mutants as misses would understate every
oracle, since no correct oracle can flag them.
"""
import numpy as np

from qcore import (MUTATIONS, oracle_exact_counts, oracle_operator,
                   oracle_state, random_circuit, simulate, wilson)

rng = np.random.default_rng(20260901)

N, DEPTH, TRIALS = 4, 14, 2000

FAMILIES = {
    "{H,X,Z,Ry,Rz,CX,CZ}": ["H", "X", "Z", "Ry", "Rz", "CX", "CZ"],
    "+ T": ["H", "X", "Z", "Ry", "Rz", "CX", "CZ", "T"],
    "+ Rx": ["H", "X", "Z", "Ry", "Rz", "CX", "CZ", "Rx"],
}

mutate = MUTATIONS["Rz convention (all gates)"]

print(f"{TRIALS} circuits per family, {N} qubits, depth {DEPTH}")
print("fault: every Rz built with the opposite sign")
print("power = detected / non-equivalent mutants\n")
print(f"{'family':<24}{'equiv':>7}{'counts':>11}{'state':>10}{'operator':>11}")
print("-" * 63)

results = {}
for name, gates in FAMILIES.items():
    live = equiv = c = s = 0
    for _ in range(TRIALS):
        ops = random_circuit(rng, N, DEPTH, gates)
        mutant = mutate(ops, rng)
        if mutant is None:
            continue
        Ua, Ub = simulate(ops, N), simulate(mutant, N)
        if not oracle_operator(Ua, Ub):
            equiv += 1
            continue
        live += 1
        c += oracle_exact_counts(Ua, Ub)
        s += oracle_state(Ua, Ub)
    results[name] = (c, live)
    eq_pct = equiv / (live + equiv) * 100 if live + equiv else 0.0
    print(f"{name:<24}{eq_pct:>6.1f}%{c/live*100:>10.1f}%"
          f"{s/live*100:>9.1f}%{100.0:>10.1f}%")

print("\nThe operator oracle is the reference: equivalence is defined by it,")
print("so its power is 1 by construction, not a measured result.\n")

k, n = results["{H,X,Z,Ry,Rz,CX,CZ}"]
p, lo, hi = wilson(k, n)
print(f"95% Wilson interval on the blind family, counts oracle:")
print(f"  {k}/{n} = {p*100:.1f}%  [{lo*100:.2f}%, {hi*100:.2f}%]")
