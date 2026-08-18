#!/usr/bin/env python3
"""The generator's RNG stream is part of the bug-zoo contract.

Every zoo entry records a generator_seed, and bug-zoo/README.md calls the
entries reproducible. That only holds if a given seed keeps producing the same
circuit, so the number of RNG draws per gate cannot change silently.

This caught a real regression: adding boundary angles introduced an
unconditional rng.random() inside random_angle, which preserved the
distribution but shifted the sequence. Seed 2358, recorded in
qiskit-16594-single-qubit-003.json, stopped regenerating its circuit. The
divergence started at the second gate.

Run: python3 tools/test_generator_stream.py
"""
import math
import random
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from generator import SPECIAL_ANGLES, random_angle, random_ops  # noqa: E402


def uniform_only(rng: random.Random, num_qubits: int, depth: int) -> list:
    """What random_ops produced before boundary angles existed.

    One draw per parameter, straight from the uniform. Kept as a literal
    reference rather than a call into the current code, so that a change to the
    current code cannot quietly redefine what we are comparing against.
    """
    from generator import GATES  # noqa: PLC0415

    eligible = [g for g in GATES if g[1] <= num_qubits]
    weights = [g[3] for g in eligible]
    ops = []
    for _ in range(depth):
        name, nq, npar, _ = rng.choices(eligible, weights=weights)[0]
        qubits = rng.sample(range(num_qubits), nq)
        params = [rng.uniform(-math.pi, math.pi) for _ in range(npar)]
        ops.append((name, qubits, params))
    return ops


def check(label: str, condition: bool) -> bool:
    print(f"  [{'OK' if condition else 'FALHA'}] {label}")
    return condition


def main() -> int:
    ok = True

    print("stream inalterado quando o vies esta desligado")
    for seed in (2358, 42, 777, 1, 999):
        for depth in (5, 10, 20, 40):
            a = uniform_only(random.Random(seed), 4, depth)
            b = random_ops(random.Random(seed), 4, depth, boundary_prob=0.0)
            if a != b:
                ok = check(f"seed={seed} depth={depth}", False)
                break
    ok = check("todas as combinacoes de seed e depth", ok) and ok

    print("o vies ainda funciona quando ligado")
    hits = sum(
        1
        for i in range(200)
        if any(
            abs(random_angle(random.Random(7 + i), boundary_prob=1.0) - s) < 1e-15
            for s in SPECIAL_ANGLES
        )
    )
    ok = check(f"boundary_prob=1.0 cai em SPECIAL_ANGLES ({hits}/200)", hits == 200) and ok

    misses = sum(
        1
        for i in range(200)
        if any(
            abs(random_angle(random.Random(7 + i), boundary_prob=0.0) - s) < 1e-15
            for s in SPECIAL_ANGLES
        )
    )
    ok = check(f"boundary_prob=0.0 nunca cai neles ({misses}/200)", misses == 0) and ok

    print("OK" if ok else "FALHOU")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
