"""Differential test of stim's own circuit transformations.

Result on 2026-09-16, stim 1.16.0: 400 random circuits, four transformations,
zero divergences. Requires numpy and stim.

That number means nothing on its own, so the second half of this file breaks the
transformations on purpose and reports what the oracle catches: 84 of 200 for a
dropped unitary, 117 of 156 for reversed CX operands. Those are raw rates rather
than powers, because these circuits are not deterministic and there is no exact
reference here to separate an equivalent mutant from a missed one. So: clean,
under an oracle with declared and partial sensitivity.

Every transformation below is documented to preserve what the circuit computes.

The oracle is the distribution of measurement outcomes over many shots, not a
single record. A first version used one shot on the assumption that a noiseless
Clifford circuit measures deterministically. That is false as soon as a
Hadamard precedes a measurement, and it produced 145 phantom divergences that
were nothing but different random draws.

inverse() is not in the list. stim refuses to invert a circuit containing resets
or measurements, correctly, and every circuit generated here has both.
"""
import numpy as np, stim

rng = np.random.default_rng(20260916)
ONE = ["H", "X", "Y", "Z", "S", "S_DAG", "SQRT_X"]
TWO = ["CX", "CZ", "SWAP", "ISWAP", "XCZ", "YCZ"]

def random_circuit(n, depth, with_repeat):
    c = stim.Circuit()
    c.append("R", list(range(n)))
    flipped = [q for q in range(n) if rng.random() < 0.4]
    if flipped: c.append("X", flipped)
    body = stim.Circuit()
    for _ in range(depth):
        if n >= 2 and rng.random() < 0.5:
            a, b = rng.choice(n, 2, replace=False)
            body.append(str(rng.choice(TWO)), [int(a), int(b)])
        else:
            body.append(str(rng.choice(ONE)), [int(rng.integers(n))])
        if rng.random() < 0.3:
            body.append("M", [int(rng.integers(n))])
    if with_repeat and len(body) > 0:
        c += body * int(rng.integers(2, 4))
    else:
        c += body
    c.append("M", list(range(n)))
    return c

SHOTS = 20_000
TOL = 1.2e-2

def record(c):
    try:
        return c.compile_sampler(seed=12345).sample(SHOTS).mean(axis=0)
    except Exception as e:
        return ("ERROR", type(e).__name__, str(e)[:60])

def differs(a, b):
    if a.shape != b.shape:
        return True
    return bool(np.abs(a - b).max() > TOL)

TRANSFORMS = {
    "flattened()":        lambda c: c.flattened(),
    "str round-trip":     lambda c: stim.Circuit(str(c)),
    "without_noise()":    lambda c: c.without_noise(),
    "copy via +":         lambda c: stim.Circuit() + c,
}

TRIALS = 400
bad = {k: [] for k in TRANSFORMS}
applied = {k: 0 for k in TRANSFORMS}

for t in range(TRIALS):
    n = int(rng.integers(2, 6))
    c = random_circuit(n, int(rng.integers(3, 10)), with_repeat=(t % 2 == 0))
    r0 = record(c)
    if isinstance(r0, tuple):
        continue
    for name, fn in TRANSFORMS.items():
        try:
            c2 = fn(c)
        except Exception as e:
            bad[name].append((str(c), f"transform raised {type(e).__name__}: {e}"))
            continue
        applied[name] += 1
        r1 = record(c2)
        if isinstance(r1, tuple):
            bad[name].append((str(c), f"transformed circuit failed to sample: {r1[1]}"))
        elif differs(r0, r1):
            d = np.abs(r0 - r1).max() if r0.shape == r1.shape else float("nan")
            bad[name].append((str(c), f"distribution moved by {d:.4f}"))

print(f"{TRIALS} circuitos aleatorios, oraculo = distribuicao de medicao, 20k shots\n")
print(f"{'transformacao':<22} {'aplicada':>9} {'divergencias':>13}")
print("-" * 48)
for name in TRANSFORMS:
    print(f"{name:<22} {applied[name]:>9} {len(bad[name]):>13}")
print("-" * 48)
for name, lst in bad.items():
    if lst:
        print(f"\n### {name}: {len(lst)} caso(s). Primeiro:")
        circ, why = lst[0]
        print("  ", why)
        print("  circuito:\n" + "\n".join("    " + l for l in circ.splitlines()[:14]))

# ---------------------------------------------------------------------------
# Does this test have any power? A clean run is worthless without knowing that.
# These transformations are wrong on purpose. Anything the oracle misses here is
# a fault class the clean run above never had the ability to see.
# ---------------------------------------------------------------------------
def mut_drop_last_gate(c):
    ins = [c[i] for i in range(len(c))]
    spots = [i for i, x in enumerate(ins) if x.name in ONE + TWO]
    if not spots: return c
    del ins[spots[-1]]
    out = stim.Circuit()
    for x in ins: out.append(x)
    return out

def mut_swap_two_qubit_order(c):
    ins = [c[i] for i in range(len(c))]
    out = stim.Circuit()
    for x in ins:
        if x.name in ("CX", "XCZ", "YCZ"):
            t = [q.value for q in x.targets_copy()]
            out.append(x.name, t[::-1])
        else:
            out.append(x)
    return out

MUTANTS = {"drop last unitary": mut_drop_last_gate,
           "reverse CX operands": mut_swap_two_qubit_order}

print()
print("controle: transformacoes QUEBRADAS de proposito")
print(f"{'mutante':<24} {'aplicado':>9} {'pego pelo oraculo':>19}")
print("-" * 56)
caught = {k: 0 for k in MUTANTS}; tried = {k: 0 for k in MUTANTS}
rng2 = np.random.default_rng(7)
for t in range(200):
    n = int(rng2.integers(2, 6))
    c = random_circuit(n, int(rng2.integers(3, 10)), with_repeat=(t % 2 == 0))
    r0 = record(c)
    if isinstance(r0, tuple): continue
    for name, fn in MUTANTS.items():
        try: c2 = fn(c.flattened())
        except Exception: continue
        if str(c2) == str(c.flattened()): continue
        tried[name] += 1
        r1 = record(c2)
        if isinstance(r1, tuple) or differs(r0, r1): caught[name] += 1
for name in MUTANTS:
    pct = f"{caught[name]}/{tried[name]}" if tried[name] else "n/a"
    print(f"{name:<24} {tried[name]:>9} {pct:>19}")
