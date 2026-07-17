#!/usr/bin/env python3
"""Classify a minimal witness as a known bug class, or as something new.

A campaign against Qiskit 2.5.0 keeps rediscovering #16594: it is reachable in
roughly 1.5% of random circuits, so a 3000-circuit run yields tens of instances
that all shrink to the same core. Reporting that as a bug count would be a lie
of the kind this project has already paid for once. The number that matters is
how many findings are NOT the known class.

The classifier computes the accumulated rotation rather than matching gate
names. A name-based rule ("two sxdg on a qubit plus a gate elsewhere") would
swallow a genuinely new bug that happens to share that shape, which is the
opposite of what this is for. Being too eager to call something known is the
dangerous direction here, so the angle condition is checked exactly.
"""

from __future__ import annotations

import math
from collections import Counter

# X-rotation family and the angle each contributes.
_X_ANGLE = {
    "sx": math.pi / 2,
    "sxdg": -math.pi / 2,
    "x": math.pi,
}

_TOL = 1e-9


def _gate_fields(g):
    if isinstance(g, dict):
        return g["gate"], g["qubits"], g.get("params", [])
    return g[0], g[1], g[2]


def _x_rotation(name, params) -> float | None:
    """Angle if this is an X-family rotation, else None."""
    if name in _X_ANGLE:
        return _X_ANGLE[name]
    if name == "rx":
        return float(params[0])
    return None


def _is_negative_half_pi_multiple(theta: float) -> bool:
    """The #16594 condition, measured: a negative multiple of pi/2 that is not
    a multiple of 2*pi. Even multiples count -- the maintainer's description
    said odd, but -pi and -3pi break too, and 2pi multiples only survive
    because emitting nothing is correct there.
    """
    if theta >= -_TOL:
        return False
    k = theta / (math.pi / 2)
    if abs(k - round(k)) > 1e-6:
        return False
    return round(k) % 4 != 0


def classify(gates) -> str:
    """Known class name, or 'unclassified' for something worth reading."""
    per: dict[int, list[tuple[str, list]]] = {}
    for g in gates:
        name, qubits, params = _gate_fields(g)
        for q in qubits:
            per.setdefault(q, []).append((name, params))

    # Qiskit #16594: one qubit carries only X-family rotations summing to a
    # negative multiple of pi/2, and at least one operation sits on another
    # qubit to push the pass onto the commutation path.
    if len(per) >= 2:
        for q, seq in per.items():
            angles = [_x_rotation(n, p) for n, p in seq]
            if not angles or any(a is None for a in angles):
                continue
            if len(angles) < 2:
                continue
            if _is_negative_half_pi_multiple(sum(angles)):
                return "qiskit-16594-commutative-cancellation"

    # A lone controlled-phase at a small angle: synthesis contract, not a bug.
    if len(gates) == 1:
        name, _, params = _gate_fields(gates[0])
        if name in ("cp", "crz") and params and abs(params[0]) < 0.1:
            return "qiskit-cp-small-angle-contract"

    return "unclassified"


def summarize(entries) -> Counter:
    return Counter(classify(e["gates"]) for e in entries)


if __name__ == "__main__":
    import json
    import pathlib
    import sys

    zoo = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "bug-zoo")
    entries = []
    for f in sorted(zoo.glob("*.json")):
        try:
            d = json.loads(f.read_text())
            if "gates" in d:
                entries.append(d)
        except Exception:  # noqa: BLE001
            continue

    counts = summarize(entries)
    print(f"{len(entries)} zoo entries")
    for name, n in counts.most_common():
        marker = "   <-- READ THESE" if name == "unclassified" else ""
        print(f"  {n:>4}  {name}{marker}")

    new = [e for e in entries if classify(e["gates"]) == "unclassified"]
    if new:
        print("\nUnclassified witnesses:")
        for e in new:
            gates = " ".join(
                f"{_gate_fields(g)[0]}{_gate_fields(g)[1]}" for g in e["gates"]
            )
            print(f"  {e['id']}: {gates}")
